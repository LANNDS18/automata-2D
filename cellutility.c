#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"

/*
 * Initialize with the fraction of filled cells equal to rho
 */

void init_cell_with_seed(int l, int seed, double rho, int *ncell, int **allcell) {
    rinit(seed);
    double r;
    *ncell = 0;
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < l; j++) {
            r = uni();

            if (r < rho) {
                allcell[i][j] = 1;
                (*ncell)++;
            } else {
                allcell[i][j] = 0;
            }
        }
    }
}

/*
 * Initialize the assigning part of allcell to cell by lx, ly and their start point's coord
 */

void init_local_cell(int lx, int ly, int *coord, int **allcell, int **cell) {
    int i, j;

    /* Copy the result from allcell based on COORD with LX and LY for local simulation */
    for (i = 1; i <= lx; i++) {
        for (j = 1; j <= ly; j++) {
            cell[i][j] = allcell[coord[0] + i - 1][coord[1] + j - 1];
        }
    }

    for (i = 0; i <= lx + 1; i++) // zero the bottom and top halos
    {
        cell[i][0] = 0;
        cell[i][ly + 1] = 0;
    }

    for (j = 0; j <= ly + 1; j++) // zero the left and right halos
    {
        cell[0][j] = 0;
        cell[lx + 1][j] = 0;
    }
}


/*
 * Display a 2d-array by lx and ly, useful for debugging
 */

void print_2d_array(int lx, int ly, int **display) {
    for (int j = ly - 1; j >= 0; j--) {
        for (int i = 0; i < lx; i++) {
            printf("%d,\t", display[i][j]);
        }
        printf("\n");
    }
}

/*
 * Print the timing and updating results
 */

void print_updating_result(double t_end, double t_start, int step, int ncell, int upper_target, int lower_target,
                           int maxstep) {
    double interval = t_end - t_start;
    double time_per_step = interval / step;

    printf("*************************************************************************************\n");
    printf("cellutility: Total computing time is %f [s]\n", interval);
    printf("cellutility: Time per step is %g [s]\n", time_per_step);
    printf("\n");
    if (ncell >= upper_target) {
        printf("cellutility: Sucesslly achieve upper target, current live cells: %d, step: %d\n", ncell, step);
    } else if (ncell <= lower_target) {
        printf("cellutility: Sucesslly achieve lower target, current live cells: %d, step: %d\n", ncell, step);
    } else {
        printf("cellutility: Fail to achieve target, exceed max steps:  %d, current live cells: %d\n", maxstep, ncell);
    }
    printf("*************************************************************************************\n");
}

/*
 * Compute number of live cells within a cell by lx and ly
 * Return int total: Number of total live cell
 */

int check_number_live_cells(int lx, int ly, int **cell) {
    int total = 0;
    for (int i = 0; i < lx; i++) {
        for (int j = 0; j < ly; j++) {
            if (cell[i][j] == 1)
                total++;
        }
    }
    printf("debug: Number of Total Live Cells: %d\n", total);
    return total;
}

/*
 * Check the count and value of the argument
 */

int check_argument(int argc, char *argv[], int *seed, int *L, double *rho) {
    switch (argc) {
        case 2:
            printf("args: L is set to default (768), rho is set to default (0.49)\n");
            *seed = atoi(argv[1]);
            *L = 768;
            *rho = 0.49;
            break;
        case 3:
            printf("args: rho is set to default (0.49)\n");
            *seed = atoi(argv[1]);
            *L = atoi(argv[2]);
            *rho = 0.49;
            break;
        case 4:
            *seed = atoi(argv[1]);
            *L = atoi(argv[2]);
            *rho = atof(argv[3]);
            break;
        default:
            printf("args: invalid\n");
            printf("Usage: automaton <int: seed (must be specified)> <int: L (default 768)> <double: rho (default 0.49)>\n");
            return 0;
    }

    if (*L < 1) {
        printf("args: invalid input, L must be greater than 0\n");
        return 0;
    }

    if (*rho <= 0 || *rho >= 1) {
        printf("args: invalid input, rho must be greater than 0 and less than 1\n");
        return 0;
    }
    return 1;
}


/*
 * print information after initializing cells and 2d grid for decomposition
 */
void print_init_cell_info(int L, double rho, int seed, int maxstep, int size, int ncell, int l_target, int u_target,
                          int *dim) {
    int LX = (L / dim[0]);
    int LY = (L / dim[1]);
    printf("-----------------------------------------------------------------------\n");

    printf("cell init info: L = %d, rho = %f, seed = %d, maxstep = %d\n",
           L, rho, seed, maxstep);
    printf("cell init info: running on %d process(es)\n", size);
    printf("cell init info: 2D decomposition grid: (%d, %d)\n", dim[0], dim[1]);

    printf("cell init info: LX= %d, LY=%d\n", LX, LY);
    if (L % dim[0] > 0)
        printf("cell init info: LX = %d, for processes allocated at (%d,*) in 2d cartesian\n",
               L - (L / dim[0]) * (dim[0] - 1), dim[0]);
    if (L % dim[1] > 0)
        printf("cell init info: LY = %d, for processes allocated at (*,%d) in 2d cartesian\n",
               L - (L / dim[1]) * (dim[1] - 1), dim[1]);

    printf("cell init info: rho = %f, live cells = %d, actual density = %f\n",
           rho, ncell, ((double) ncell) / ((double) L * L));

    printf("cell init info: lower target number of cells: %d\n", l_target);
    printf("cell init info: upper target number of cells: %d\n", u_target);

    if (u_target > L * L)
        printf("cell init info: upper target equals to %d which is unachievable\n", u_target);
    if (l_target < 0)
        printf("cell init info: upper target equals to %d which is unachievable\n", l_target);

    printf("-----------------------------------------------------------------------\n");
}
