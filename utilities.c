#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"

void initial_array_with_0(int lx, int ly, int **cell)
{
    for (int i = 0; i < lx; i++)
    {
        for (int j = 0; j < ly; j++)
        {
            cell[i][j] = 0;
        }
    }
}

void print_array(int lx, int ly, int **display)
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