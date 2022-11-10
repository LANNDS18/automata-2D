/*
 *  Main header file for percolation exercise.
 */

/*
 *  System size L
 */

#define L 500

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
void initial_array_with_0(int l, int cell[l][l]);
void print_array(int l, int display[l][l]);
void check_number_live_cells(int l, int cell[l][l]);