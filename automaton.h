#include <mpi.h>

/*
 *  System size L
 */

#define L 768

#define NPROC 4

#define TRUE 1
#define FALSE 0

struct adjacent_process
{
    int left;
    int right;
    int up;
    int down;
};
/*
 * Cell updating functions implemented with MPI
 */

void create_2d_cart_and_assign_coord(int rank, MPI_Comm comm, MPI_Comm *cart, int *LX, int *LY, int *COORD);
struct adjacent_process get_adjacent_processes(MPI_Comm cart);
void halo_swap_2d_mpi(int lx, int ly, int **cell, struct adjacent_process p_adj, MPI_Comm cart, MPI_Datatype VERTICAL_HALO_TYPE, MPI_Request *request);
int update_live_cell_mpi(int lx, int ly, int **neigh, int **cell, MPI_Request *request);


/*
 *  Visualisation
 */

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
int check_number_live_cells(int lx, int ly, int **cell);

/*
 * Dynamic Array Allocation
 */

void *arralloc(size_t size, int ndim, ...);