/*
 *  Main header file for percolation exercise.
 */

/*
 *  System size L
 */

#define L 768

#define NPROC 4

#define TRUE 1
#define FALSE 0

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
 * utilities function for initialization and debugging
 */
void initial_array_with_0(int l, int **cell);
void print_array(int l, int **display);
void check_number_live_cells(int l, int **cell);

/*
 * Dynamic Array Allocation 
 */

void *arralloc(size_t size, int ndim, ...);