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
    for (int i = 0; i < lx; i++)
    {
        for (int j = 0; j < ly; j++)
        {
            printf("%d,", display[i][j]);
        }
        printf("\n");
    }
}


/*
 * Compute number of live cells within a cell by lx and ly
 */

void check_number_live_cells(int lx, int ly, int **cell)
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
}