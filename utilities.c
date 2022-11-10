#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"

void initial_array_with_0(int l, int **cell)
{
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            cell[i][j] = 0;
        }
    }
}

void print_array(int l, int **display)
{
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            printf("%d,", display[i][j]);
        }
        printf("\n");
    }
}

void check_number_live_cells(int l, int **cell)
{
    int total = 0;
    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            if (cell[i][j] == 1)
                total++;
        }
    }
    printf("Number of Total Live Cells: %d\n", total);
}