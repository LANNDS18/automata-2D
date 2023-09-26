#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

int main(int argc, char *argv[]) {

    /*
     *  Variables that define the simulation
     */

    int seed, L, LX, LY; // LX, LY: The length assigned cell for each dim
    double rho;
    int lower_target, upper_target;

    /*
     *  Local variables
     */
    int ncell, localncell, maxstep, printfreq;
    int step = 0;
    int cell_coord[2];               // the start point of the assigned part of cell for each dim
    struct adjacent_process p_neigh; // Ranks of neighbour process in 2D

    /*
     *  MPI common world variables
     */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart;          // 2d Cart topology comm world
    MPI_Request request[4]; // Store Isend requests
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int arg_check_bit;

    if (rank == 0) {
        arg_check_bit = check_argument(argc, argv, &seed, &L, &rho);
    }

    MPI_Bcast(&arg_check_bit, 1, MPI_INT, 0, comm);

    if (arg_check_bit == 0) {
        MPI_Finalize();
        return 1;
    }

    // Brocast valur of L andto all processes
    MPI_Bcast(&L, 1, MPI_INT, 0, comm);

    /*
     *  Update for a large number of steps to prevent execute without stopping
     *  periodically report progress
     */

    maxstep = 1000000;
    printfreq = 500;

    /*
     *  Additional array WITHOUT halos for initialisation and IO.
     */

    int **allcell = arralloc(sizeof(int), 2, L, L);

    /*
     *Initalize the Dim and Create 2d Cart
     */

    int dim[2] = {0, 0};
    create_2d_cart_mpi(size, comm, dim, &cart);
    allocate_cells_mpi(L, rank, dim, cart, &LX, &LY, cell_coord);

    /*
     * Define a vector derived datatype for send and recv the halo from
     */

    MPI_Datatype VERTICAL_HALO_TYPE;
    MPI_Type_vector(LX, 1, LY + 2, MPI_INT, &VERTICAL_HALO_TYPE);
    MPI_Type_commit(&VERTICAL_HALO_TYPE);

    /*
     *  Define the main arrays for the simulation based on LX, LY
     */

    int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);
    int **neigh = arralloc(sizeof(int), 2, LX + 2, LY + 2);

    /* Get the neighbour processes */
    p_neigh = get_adjacent_processes_mpi(cart);

    if (rank == 0) {
        /*
         *  Set the randum number seed and initialise the generator
         */

        /*  Initialise with the fraction of filled cells equal to rho */
        init_cell_with_seed(L, seed, rho, &ncell, allcell);

        /*  Compute Lower and Upper Bound of Target value */
        lower_target = (int) ((double) ncell * 2 / 3);
        upper_target = (int) ((double) ncell * 3 / 2);

        print_init_cell_info(L, rho, seed, maxstep, size, ncell, lower_target, upper_target, dim);
    }

    MPI_Bcast(&allcell[0][0], L * L, MPI_INT, 0, comm);
    MPI_Bcast(&lower_target, 1, MPI_INT, 0, comm);
    MPI_Bcast(&upper_target, 1, MPI_INT, 0, comm);

    init_local_cell(LX, LY, cell_coord, allcell, cell);

    // Start the timer
    MPI_Barrier(cart);
    double t_start = MPI_Wtime();

    while (step <= maxstep) {
        step++;

        halo_swap_2d_mpi(LX, LY, cell, p_neigh, cart, VERTICAL_HALO_TYPE, request);

        /*
         * Update live cell by counting neighbour
         */

        localncell = update_live_cell_mpi(LX, LY, neigh, cell, request);

        /*
         *  Compute global number of changes on all processes for checking whether reaching target value
         */

        MPI_Allreduce(&localncell, &ncell, 1, MPI_INT, MPI_SUM, cart);

        /*
         *  Report progress by printfreq
         */

        if (step % printfreq == 0) {
            if (rank == 0) {
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
    /*
     * Output the timing result and whether and which targets it achieve
     */

    if (rank == 0) {
        print_updating_result(t_end, t_start, step, ncell, upper_target, lower_target, maxstep);
    }

    collect_allcells_mpi(L, LX, LY, cell_coord, cell, cart, allcell);

    /*
     *  Write the cells to the file "cell.pbm" from rank 0 with dynamic allocation
     */

    if (rank == 0) {
        cellwritedynamic("cell.pbm", allcell, L);
    }

    /* Free the dynamic assigned memory and Finalize */
    free(neigh);
    free(cell);
    free(allcell);

    MPI_Finalize();

    return 0;
}
