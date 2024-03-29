#include <mpi.h>

#define TRUE 1
#define FALSE 0

struct adjacent_process {
    int left;
    int right;
    int up;
    int down;
};

/*
 * Cell updating functions implemented with MPI
 */

void create_2d_cart_mpi(int size, MPI_Comm comm, int *dim, MPI_Comm *cart);

void allocate_cells_mpi(int L, int rank, int *dim, MPI_Comm cart, int *LX, int *LY, int *COORD);

struct adjacent_process get_adjacent_processes_mpi(MPI_Comm cart);

void halo_swap_2d_mpi(int lx, int ly, int **cell, struct adjacent_process p_adj, MPI_Comm cart,
                      MPI_Datatype VERTICAL_HALO_TYPE, MPI_Request *request);

int update_live_cell_mpi(int lx, int ly, int **neigh, int **cell, MPI_Request *request);

void collect_allcells_mpi(int L, int LX, int LY, int *COORD, int **cell, MPI_Comm cart, int **allcell);


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

void init_local_cell(int lx, int ly, int *coord, int **allcell, int **cell);

void print_init_cell_info(int L, double rho, int seed, int maxstep, int size, int ncell, int l_target, int u_target,
                          int *dim);

void print_updating_result(double t_end, double t_start, int step, int ncell, int upper_target, int lower_target,
                           int maxstep);

int check_argument(int argc, char *argv[], int *seed, int *L, double *rho);

/* Debugging */
void print_2d_array(int lx, int ly, int **display);

int check_number_live_cells(int lx, int ly, int **cell);

/*
 * Dynamic Array Allocation
 */

void *arralloc(size_t size, int ndim, ...);