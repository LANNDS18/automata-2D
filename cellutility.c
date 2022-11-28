#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"

/*
 * Initialize with the fraction of filled cells equal to rho
 */

void init_cell_with_seed(int l, int seed, double rho, int *ncell, int **allcell)
{
    rinit(seed);
    double r;
    *ncell = 0;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            r = uni();

            if (r < rho)
            {
                allcell[i][j] = 1;
                (*ncell)++;
            }
            else
            {
                allcell[i][j] = 0;
            }
        }
    }
}

/*
 * Initialize the assigning part of allcell to cell by lx, ly and their start point's coord
 */

void init_local_cell(int lx, int ly, int *coord, int **allcell, int **cell)
{
    int i, j;

    /* Copy the result from allcell based on COORD with LX and LY for local simulation */
    for (i = 1; i <= lx; i++)
    {
        for (j = 1; j <= ly; j++)
        {
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
 * Initialize cell by filling 0, useful for debugging and initalize temp variables
 */

void init_cell_with_0(int lx, int ly, int **cell)
{
    for (int i = 0; i < lx; i++)
    {
        for (int j = 0; j < ly; j++)
        {
            cell[i][j] = 0;
        }
    }
}

/*
 * Display a 2d-array by lx and ly
 */

void print_2d_array(int lx, int ly, int **display)
{
    for (int j = ly - 1; j >= 0; j--)
    {
        for (int i = 0; i < lx; i++)
        {
            printf("%d,", display[i][j]);
        }
        printf("\n");
    }
}

/*
 * Print the timing and updating results
 */

void print_updating_result(double t_end, double t_start, int step, int ncell, int upper_target, int lower_target, int maxstep)
{
    double interval = t_end - t_start;
    double time_per_step = interval / step;

    printf("*************************************************************************************\n");
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
        printf("Fail to achieve target, exceed max steps:  %d, current live cells: %d\n", maxstep, ncell);
    }
    printf("*************************************************************************************\n");
}

/*
 * Compute number of live cells within a cell by lx and ly
 * Return int total: Number of total live cell
 */

int check_number_live_cells(int lx, int ly, int **cell)
{
    int total = 0;
    for (int i = 0; i < lx; i++)
    {
        for (int j = 0; j < ly; j++)
        {
            if (cell[i][j] == 1)
                total++;
        }
    }
    printf("Number of Total Live Cells: %d\n", total);
    return total;
}

/*
 * Check the count and value of the argument
 */

int check_argument(int argc, char *argv[], int *seed, int *L, double *rho)
{
    switch (argc)
    {
    case 2:
        printf("L is set to default (768), rho is set to default (0.49)\n");
        *seed = atoi(argv[1]);
        *L = 768;
        *rho = 0.49;
        break;
    case 3:
        printf("rho is set to default (0.49)\n");
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
        printf("Usage: automaton <int: seed (must be specified)> <int: L (default 768)> <double: rho (default 0.49)>\n");
        return 0;
    }

    if (*L <= 0)
    {
        printf("L must be greater than 0\n");
        return 0;
    }

    if (*rho <= 0 || *rho >= 1)
    {
        printf("rho must be greater than 0 and less than 1\n");
        return 0;
    }
    return 1;
}