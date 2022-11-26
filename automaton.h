#include <mpi.h>

/*
 *  System size L
 */

#define L 768

#define NPROC 4

#define TRUE 1
#define FALSE 0


/*
 * create 2d simulation world and return assigned lx, ly coordinate to each process
 */

void create_2d_cart_and_assign_coord(int rank, MPI_Comm comm, MPI_Comm *cart, int *LX, int *LY, int *COORD);


/*
 *  Visualisation
 */

void cellwrite(char *cellfile, int cell[L][L]);
void cellwritedynamic(char *cellfile, int **cell, int l);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);


/*
 * utilities function for initialization, copy, and couting live cells
 */

void init_cell_with_seed(int l, int seed, double rho, int *ncell, int **allcell);
void init_cell_with_0(int lx, int ly, int **cell);
void init_local_cell(int lx, int ly, int *coord, int **allcell, int **cell);

void print_2d_array(int lx, int ly, int **display);
void check_number_live_cells(int lx, int ly, int **cell);

/*
 * Dynamic Array Allocation
 */

void *arralloc(size_t size, int ndim, ...);