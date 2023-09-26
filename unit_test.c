#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

void test_halo_swap() {

    int L, LX, LY; // LX, LY: The length assigned cell for each dim
    L = 4;

    int cell_coord[2];               // the start point of the assigned part of cell for each dim
    struct adjacent_process p_neigh; // Ranks of neighbour process in 2D

    /*
     *  MPI common world variables
     */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart;          // 2d Cart topology comm world
    MPI_Request request[4]; // Store Isend requests
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        printf("=============================================== \n");
        printf("Start to run Unit Test for Halo Swapping \n");
        printf("------------------------------------------------\n");
    }

    MPI_Bcast(&L, 1, MPI_INT, 0, comm);
    int **allcell = arralloc(sizeof(int), 2, L, L);

    int dim[2] = {0, 0};
    create_2d_cart_mpi(size, comm, dim, &cart);
    allocate_cells_mpi(L, rank, dim, cart, &LX, &LY, cell_coord);

    MPI_Datatype VERTICAL_HALO_TYPE;
    MPI_Type_vector(LX, 1, LY + 2, MPI_INT, &VERTICAL_HALO_TYPE);
    MPI_Type_commit(&VERTICAL_HALO_TYPE);

    int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);

    p_neigh = get_adjacent_processes_mpi(cart);

    int loop = 0;

    if (rank == 0) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                allcell[i][j] = loop;
                loop++;
            }
        }
    }

    MPI_Bcast(&allcell[0][0], L * L, MPI_INT, 0, comm);

    init_local_cell(LX, LY, cell_coord, allcell, cell);

    if (rank == 0) {
        printf("all cell:\n");
        print_2d_array(L, L, allcell);

        printf("rank 0 before swapping: \n");
        print_2d_array(LX + 2, LY + 2, cell);
    }

    halo_swap_2d_mpi(LX, LY, cell, p_neigh, cart, VERTICAL_HALO_TYPE, request);
    MPI_Status status;
    MPI_Wait(&request[0], &status);
    MPI_Wait(&request[1], &status);
    MPI_Wait(&request[2], &status);
    MPI_Wait(&request[3], &status);

    if (rank == 0) {

        printf("rank 0 after swapping: \n");
        print_2d_array(LX + 2, LY + 2, cell);
        printf("\n");
        printf("------------------------------------------------\n");
        printf("Finish Unit Test for Halo Swapping \n");
        printf("=============================================== \n");
    }

    free(cell);
    free(allcell);
}

void test_2d_grid() {
    int L, LX, LY; // LX, LY: The length assigned cell for each dim
    L = 7;

    int cell_coord[2]; // the start point of the assigned part of cell for each dim

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart; // 2d Cart topology comm world
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        printf("=============================================== \n");
        printf("Start to run Unit Test for cartesian topology \n");
        printf("------------------------------------------------\n");
    }

    MPI_Bcast(&L, 1, MPI_INT, 0, comm);

    int dim[2] = {0, 0};
    create_2d_cart_mpi(size, comm, dim, &cart);
    allocate_cells_mpi(L, rank, dim, cart, &LX, &LY, cell_coord);
    struct adjacent_process p_neigh = get_adjacent_processes_mpi(cart);
    MPI_Barrier(comm);
    if (rank == 0) {
        printf("rank %d, left: %d, right: %d, up: %d, down: %d, LX:%d, LY:%d\n", rank, p_neigh.left, p_neigh.right,
               p_neigh.up, p_neigh.down, LX, LY);
        printf("------------------------------------------------\n");
        printf("Finish Unit Test for cartesian topology \n");
        printf("=============================================== \n");
    }
}

void test_results_collection() {
    int L, LX, LY; // LX, LY: The length assigned cell for each dim
    L = 7;

    int cell_coord[2]; // the start point of the assigned part of cell for each dim

    /*
     *  MPI common world variables
     */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart; // 2d Cart topology comm world
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        printf("=============================================== \n");
        printf("Start to run Unit Test for Result Collection  \n");
        printf("------------------------------------------------\n");
    }

    MPI_Bcast(&L, 1, MPI_INT, 0, comm);

    int dim[2] = {0, 0};
    create_2d_cart_mpi(size, comm, dim, &cart);
    allocate_cells_mpi(L, rank, dim, cart, &LX, &LY, cell_coord);

    int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);
    int **allcell = arralloc(sizeof(int), 2, L, L);

    for (int i = 1; i <= LX + 1; i++) {
        for (int j = 1; j <= LY + 1; j++) {
            cell[i][j] = rank;
        }
    }

    collect_allcells_mpi(L, LX, LY, cell_coord, cell, cart, allcell);

    if (rank == 0) {
        printf("all cell:\n");
        print_2d_array(L, L, allcell);
        printf("------------------------------------------------\n");
        printf("Finish Unit Test for Result Collection \n");
        printf("=============================================== \n");
    }

    free(allcell);
    free(cell);
}

void test_unbalanced_distribution() {
    int L, LX, LY; // LX, LY: The length assigned cell for each dim
    L = 7;

    int cell_coord[2]; // the start point of the assigned part of cell for each dim

    /*
     *  MPI common world variables
     */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm cart; // 2d Cart topology comm world
    int size, rank;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        printf("=============================================== \n");
        printf("Start to run Unit Test for Unbalanced Distribution  \n");
        printf("------------------------------------------------\n");
    }

    MPI_Bcast(&L, 1, MPI_INT, 0, comm);

    int dim[2] = {0, 0};
    create_2d_cart_mpi(size, comm, dim, &cart);
    allocate_cells_mpi(L, rank, dim, cart, &LX, &LY, cell_coord);

    int **cell = arralloc(sizeof(int), 2, LX + 2, LY + 2);
    int **allcell = arralloc(sizeof(int), 2, L, L);

    for (int i = 1; i <= LX + 1; i++) {
        for (int j = 1; j <= LY + 1; j++) {
            cell[i][j] = rank;
        }
    }

    collect_allcells_mpi(L, LX, LY, cell_coord, cell, cart, allcell);

    if (rank == 0) {
        printf("all cell:\n");
        print_2d_array(L, L, allcell);
        printf("------------------------------------------------\n");
        printf("Finish Unit Test for Unbalanced Distribution \n");
        printf("=============================================== \n");
    }
    free(allcell);
    free(cell);
}

int main() {
    MPI_Init(NULL, NULL);
    test_2d_grid();
    test_halo_swap();
    test_results_collection();
    test_unbalanced_distribution();
    MPI_Finalize();
    return 0;
}
