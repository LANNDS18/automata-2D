#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

/*
 * 2D decomposition parallel program to simulate a 2D cellular automaton
 */

int main(int argc, char *argv[])
{

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  int **allcell = arralloc(sizeof(int), 2, L, L);
  int **tempcell = arralloc(sizeof(int), 2, L, L);

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */
  int lower_target, upper_target;
  int i, j, ncell, localncell, maxstep, printfreq;
  int step = 0;
  double r;

  /*
   *  MPI common world variables
   */

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;
  MPI_Request request[4];

  int size, rank;
  int tag = 1;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (NPROC != size)
  {
    if (rank == 0)
    {
      printf("automaton: ERROR, NPROC = %d running on %d process(es)\n",
             NPROC, size);
    }

    MPI_Finalize();
    return 1;
  }

  if (argc != 2)
  {
    if (rank == 0)
    {
      printf("Usage: automaton <seed>\n");
    }

    MPI_Finalize();
    return 1;
  }

  /*
   *  Update for a large number of steps to prevent execute without stopping
   *  periodically report progress
   */

  maxstep = 1000000;
  printfreq = 500;

  if (rank == 0)
  {
    printf("automaton: running on %d process(es)\n", size);

    /*
     *  Set the cell density rho (between 0 and 1)
     */

    rho = 0.49;

    /*
     *  Set the randum number seed and initialise the generator
     */

    seed = atoi(argv[1]);

    printf("automaton: L = %d, rho = %f, seed = %d, maxstep = %d\n",
           L, rho, seed, maxstep);

    rinit(seed);

    /*
     *  Initialise with the fraction of filled cells equal to rho
     */

    ncell = 0;

    for (i = 0; i < L; i++)
    {
      for (j = 0; j < L; j++)
      {
        r = uni();

        if (r < rho)
        {
          allcell[i][j] = 1;
          ncell++;
        }
        else
        {
          allcell[i][j] = 0;
        }
      }
    }
    lower_target = (int)((double)ncell * 2 / 3);
    upper_target = (int)((double)ncell * 3 / 2);

    printf("automaton: rho = %f, live cells = %d, actual density = %f\n",
           rho, ncell, ((double)ncell) / ((double)L * L));
    printf("lower target number of cells: %d\n", lower_target);
    printf("upper target number of cells: %d\n", upper_target);
  }

  // Deistribute the allcell to all the processess
  MPI_Bcast(&allcell[0][0], L * L, MPI_INT, 0, comm);
  MPI_Bcast(&lower_target, 1, MPI_INT, 0, comm);
  MPI_Bcast(&upper_target, 1, MPI_INT, 0, comm);

  initial_array_with_0(L, L, tempcell);

  // Initalize the Dim and Cart
  int coord[2];
  int dim[2] = {0, 0};
  int up, down, left, right;
  int period[2] = {FALSE, TRUE};
  int dim_size = 2;
  MPI_Comm cart;

  /* Create a 2D (n * m) Cartersian Topology where n*m = NPROC */
  MPI_Dims_create(NPROC, 2, dim);
  MPI_Cart_create(comm, dim_size, dim, period, FALSE, &cart);
  MPI_Cart_coords(cart, rank, dim_size, coord);

  int LX = L / dim[0];
  int LY = L / dim[1];

  int X_COORD = L / dim[0] * coord[0];
  int Y_COORD = L / dim[1] * coord[1];

  int IS_END_X = coord[0] == (dim[0] - 1);
  int IS_END_Y = coord[1] == (dim[1] - 1);

  if (IS_END_X)
    LX = L - LX * (dim[0] - 1);
  if (IS_END_Y)
    LY = L - LY * (dim[1] - 1);

  printf("L = %d, LY= %d, LX=%d, Rank=%d\n", L, LY, LX, rank);

  // printf("%d, %d, coord: (%d, %d), LX, LY,", IS_END_X, IS_END_Y, coord[0], coord[1]);

  MPI_Cart_shift(cart, 0, 1, &left, &right);
  MPI_Cart_shift(cart, 1, 1, &down, &up);

  MPI_Datatype VERTICAL_HALO_TYPE;
  MPI_Type_vector(LX, 1, LY + 2, MPI_INT, &VERTICAL_HALO_TYPE);
  MPI_Type_commit(&VERTICAL_HALO_TYPE);


  /*
   *  Define the main arrays for the simulation
   */
  int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);
  int **neigh = arralloc(sizeof(int), 2, LX + 2, LY + 2);

  int **smallcell = arralloc(sizeof(int), 2, LX, LY);

  for (i = 0; i < LX; i++)
  {
    for (j = 0; j < LY; j++)
    {
      smallcell[i][j] = allcell[X_COORD + i][Y_COORD + j];
    }
  }

  for (i = 1; i <= LX; i++)
  {
    for (j = 1; j <= LY; j++)
    {
      cell[i][j] = smallcell[i - 1][j - 1];
    }
  }

  for (i = 0; i <= LX + 1; i++) // zero the bottom and top halos
  {
    cell[i][0] = 0;
    cell[i][LY + 1] = 0;
  }

  for (j = 0; j <= LY + 1; j++) // zero the left and right halos
  {
    cell[0][j] = 0;
    cell[LX + 1][j] = 0;
  }

  // Start the timer
  MPI_Barrier(cart);
  double t_start = MPI_Wtime();

  while (step <= maxstep)
  {
    step++;

    MPI_Issend(&cell[LX][1], LY, MPI_INT, right, tag, cart, &request[0]);
    MPI_Issend(&cell[1][1], LY, MPI_INT, left, tag, cart, &request[1]);

    MPI_Recv(&cell[0][1], LY, MPI_INT, left, tag, cart, &status);
    MPI_Recv(&cell[LX + 1][1], LY, MPI_INT, right, tag, cart, &status);

    /* Debugging graph of cart topology for four processes
      -------
      | 1 3 |
      | 0 2 |
      -------
    */

    /* Non Blocking Send and Receive function to Upper and Lower Processes */
    MPI_Issend(&cell[1][LY], 1, VERTICAL_HALO_TYPE, up, tag, cart, &request[2]);
    MPI_Issend(&cell[1][1], 1, VERTICAL_HALO_TYPE, down, tag, cart, &request[3]);

    MPI_Recv(&cell[1][0], 1, VERTICAL_HALO_TYPE, down, tag, cart, &status);
    MPI_Recv(&cell[1][LY + 1], 1, VERTICAL_HALO_TYPE, up, tag, cart, &status);

    /* Debuging


    if (rank == 1)
    {
      print_array(LX+2, LY+2, cell);
      printf("xcoor, ycoo= (%d, %d)\n", X_COORD, Y_COORD);
      printf("left: %d, right: %d, up: %d, down: %d for rank %d step:%d\n", left, right, up, down, rank, step);
    }

    if (rank == 1)
    {
      print_array(L+2, cell);
      printf("xcoor, ycoo= (%d, %d)\n", x_coor, y_coor);
      printf("proc:%d, xcoor:%d, ycoor:%d\n",rank, x_coor, y_coor);
      printf("left: %d, right: %d, up: %d, down: %d for rank %d step:%d\n", left, right, up, down, rank, step);
    }
    */

    for (i = 1; i <= LX; i++)
    {
      for (j = 1; j <= LY; j++)
      {
        // Compute the lives cell for each assigned cell
        neigh[i][j] = cell[i][j] + cell[i][j - 1] + cell[i][j + 1] + cell[i - 1][j] + cell[i + 1][j];
      }
    }

    MPI_Wait(&request[0], &status);
    MPI_Wait(&request[1], &status);
    MPI_Wait(&request[2], &status);
    MPI_Wait(&request[3], &status);

    localncell = 0;

    for (i = 1; i <= LX; i++)
    {
      for (j = 1; j <= LY; j++)
      {

        if (neigh[i][j] == 2 || neigh[i][j] == 4 || neigh[i][j] == 5)
        {
          cell[i][j] = 1;
          localncell++;
        }
        else
        {
          cell[i][j] = 0;
        }
      }
    }

    /*
     *  Compute global number of changes on all processes for checking whether reaching target value
     */
    MPI_Allreduce(&localncell, &ncell, 1, MPI_INT, MPI_SUM, cart);

    /*
     *  Report progress every now and then
     */

    if (step % printfreq == 0)
    {
      if (rank == 0)
      {
        printf("automaton: number of live cells on step %d is %d\n",
               step, ncell);
      }
    }

    /*
     * Check if it reach the stop condition, if it is True, break the loop.
     */

    if (ncell >= upper_target || ncell <= lower_target)
    {
      if (rank == 0 && step % printfreq != 0)
      {
        printf("automaton: number of live cells on step %d is %d\n",
               step, ncell);
      }
      break;
    }
  }

  /* End of the loop stop the timing */
  MPI_Barrier(cart);
  double t_end = MPI_Wtime();
  double interval = t_end - t_start;
  double time_per_step = interval / step;

  if (rank == 0)
  {
    printf("Total computing time is %f [s]\n", interval);
    printf("Time per step is %g [s]\n", time_per_step);
  }

  /*
   *  Copy the centre of cell, excluding the halos, into allcell
   */

  for (i = 0; i < LX; i++)
  {
    for (j = 0; j < LY; j++)
    {
      tempcell[X_COORD + i][Y_COORD + j] = cell[i + 1][j + 1];
    }
  }

  /* Collect the result by reduction */
  MPI_Reduce(&tempcell[0][0], &allcell[0][0], L * L, MPI_INT, MPI_SUM, 0, cart);

  /*
   *  Write the cells to the file "cell.pbm" from rank 0 with dynamic allocation
   */

  if (rank == 0)
  {
    cellwritedynamic("cell.pbm", allcell, L);
  }

  /* Free the dynamic assigned Memory */
  free(neigh);
  free(cell);
  free(tempcell);
  free(allcell);

  /*
   * Finalise MPI before finishing
   */

  MPI_Finalize();

  return 0;
}
