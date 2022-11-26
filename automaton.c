#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

int main(int argc, char *argv[])
{

  /*
   *  Additional array WITHOUT halos for initialisation and IO.
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

    /*  Initialise with the fraction of filled cells equal to rho */
    init_cell_with_seed(L, seed, rho, &ncell, allcell);

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

  MPI_Comm cart;

  /* store the start point of the assigned part of cell for each dim */
  int COORD[2];

  /* The length assigned cell for each dim */
  int LX, LY;

  /* Rank of neighbour process */
  int up, down, left, right;

  create_2d_cart_and_assign_coord(rank, comm, &cart, &LX, &LY, &COORD[0]);
  printf("L = %d, LY= %d, LX=%d, Rank=%d\n", L, LY, LX, rank);

  /* Get the neighbour process by cart shift */
  MPI_Cart_shift(cart, 0, 1, &left, &right);
  MPI_Cart_shift(cart, 1, 1, &down, &up);

  MPI_Datatype VERTICAL_HALO_TYPE;
  MPI_Type_vector(LX, 1, LY + 2, MPI_INT, &VERTICAL_HALO_TYPE);
  MPI_Type_commit(&VERTICAL_HALO_TYPE);

  /*
   *  Define the main arrays for the simulation based on LX, LY
   */
  int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);
  int **neigh = arralloc(sizeof(int), 2, LX + 2, LY + 2);
  init_local_cell(LX, LY, COORD, allcell, cell);

  // Start the timer
  MPI_Barrier(cart);
  double t_start = MPI_Wtime();

  while (step <= maxstep)
  {
    step++;

    MPI_Isend(&cell[LX][1], LY, MPI_INT, right, tag, cart, &request[0]);
    MPI_Isend(&cell[1][1], LY, MPI_INT, left, tag, cart, &request[1]);

    MPI_Recv(&cell[0][1], LY, MPI_INT, left, tag, cart, &status);
    MPI_Recv(&cell[LX + 1][1], LY, MPI_INT, right, tag, cart, &status);

    /* Non Blocking Send and Receive function to Upper and Lower Processes */
    MPI_Isend(&cell[1][LY], 1, VERTICAL_HALO_TYPE, up, tag, cart, &request[2]);
    MPI_Isend(&cell[1][1], 1, VERTICAL_HALO_TYPE, down, tag, cart, &request[3]);

    MPI_Recv(&cell[1][0], 1, VERTICAL_HALO_TYPE, down, tag, cart, &status);
    MPI_Recv(&cell[1][LY + 1], 1, VERTICAL_HALO_TYPE, up, tag, cart, &status);

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

    /* Update live cell by counting neighbour */
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
     *  Report progress by printfreq
     */

    if (step % printfreq == 0)
    {
      if (rank == 0)
      {
        printf("automaton: number of live cells on step %d is %d\n", step, ncell);
      }
    }

    /*
     * Check if it reach the stop condition, if it is True, break the loop.
     */

    if (ncell >= upper_target || ncell <= lower_target)
      break;
  }

  /*
   *Stop the timing at the end of loop
   */

  MPI_Barrier(cart);
  double t_end = MPI_Wtime();
  double interval = t_end - t_start;
  double time_per_step = interval / step;

  /*
   * Output the timing result and whether and which targets it achieve
   */

  if (rank == 0)
  {
    printf("\n");
    printf("Total computing time is %f [s]\n", interval);
    printf("Time per step is %g [s]\n", time_per_step);
    printf("\n");
    if (ncell >= upper_target)
    {
      printf("Sucesslly achieve upper target, current live cells: %d, step: %d\n", ncell, step);
    }
    else if (ncell <= lower_target)
    {
      printf("Sucesslly achieve lower target, current live cells: %d, step: %d\n", ncell, step);
    }
    else
    {
      printf("Fail to achieve target, exceed max steps:  %d, current live cells: %d", maxstep, ncell);
    }
    printf("\n");
  }

  /* Initialize tempcell with all 0 to avoid wrong reduction */
  init_cell_with_0(L, L, tempcell);

  /*
   *  Copy the centre of cell, excluding the halos, into tempcell
   */

  for (i = 0; i < LX; i++)
  {
    for (j = 0; j < LY; j++)
    {
      tempcell[COORD[0] + i][COORD[1] + j] = cell[i + 1][j + 1];
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

  /* Free the dynamic assigned memory and Finalize */
  free(neigh);
  free(cell);
  free(tempcell);
  free(allcell);

  MPI_Finalize();

  return 0;
}

/*
 * 2D decomposition parallel program to simulate a 2D cellular automaton
 */

void create_2d_cart_and_assign_coord(int rank, MPI_Comm comm, MPI_Comm *cart, int *LX, int *LY, int *COORD)
{
  /* Initalize the Dim and Cart */
  int dim[2] = {0, 0};
  /* Set periodic in second  to TRUE */
  int period[2] = {FALSE, TRUE};
  int ndim = 2;

  /* Create a 2D (n * m) Cartersian Topology where n*m = NPROC */
  MPI_Dims_create(NPROC, 2, dim);
  MPI_Cart_create(comm, ndim, dim, period, FALSE, &(*cart));

  /* Coordinate of the current process in Cart Topology */
  int p_coord[2];
  MPI_Cart_coords(*cart, rank, ndim, p_coord);

  *LX = L / dim[0];
  *LY = L / dim[1];

  COORD[0] = L / dim[0] * p_coord[0];
  COORD[1] = L / dim[1] * p_coord[1];

  int IS_END_X = p_coord[0] == (dim[0] - 1);
  int IS_END_Y = p_coord[1] == (dim[1] - 1);

  if (IS_END_X)
    *LX = L - *LX * (dim[0] - 1);
  if (IS_END_Y)
    *LY = L - *LY * (dim[1] - 1);
}
