#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

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

void halo_swap_mpi(int LX, int LY, int **cell, int right, int left, int up, int down, MPI_Comm cart, MPI_Datatype VERTICAL_HALO_TYPE, MPI_Request *request, MPI_Status status)
{
    int tag = 1;
    /* Non Blocking Send */
    MPI_Isend(&cell[LX][1], LY, MPI_INT, right, tag, cart, &request[0]);
    MPI_Recv(&cell[0][1], LY, MPI_INT, left, tag, cart, &status);

    MPI_Isend(&cell[1][1], LY, MPI_INT, left, tag, cart, &request[1]);
    MPI_Recv(&cell[LX + 1][1], LY, MPI_INT, right, tag, cart, &status);

    MPI_Isend(&cell[1][LY], 1, VERTICAL_HALO_TYPE, up, tag, cart, &request[2]);
    MPI_Recv(&cell[1][0], 1, VERTICAL_HALO_TYPE, down, tag, cart, &status);

    MPI_Isend(&cell[1][1], 1, VERTICAL_HALO_TYPE, down, tag, cart, &request[3]);
    MPI_Recv(&cell[1][LY + 1], 1, VERTICAL_HALO_TYPE, up, tag, cart, &status);
}

int update_live_cell_mpi(int lx, int ly, int **neigh, int **cell, MPI_Request *request, MPI_Status status)
{
    int localncell = 0;
    for (int i = 1; i <= lx; i++)
    {
        for (int j = 1; j <= ly; j++)
        {
            // Compute the lives cell for each assigned cell
            neigh[i][j] = cell[i][j] + cell[i][j - 1] + cell[i][j + 1] + cell[i - 1][j] + cell[i + 1][j];
        }
    }

    MPI_Wait(&request[0], &status);
    MPI_Wait(&request[1], &status);
    MPI_Wait(&request[2], &status);
    MPI_Wait(&request[3], &status);

    for (int i = 1; i <= lx; i++)
    {
        for (int j = 1; j <= ly; j++)
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

    return localncell;
}