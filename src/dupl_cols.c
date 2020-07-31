/* $ID: dupl_cols.c, last updated 2020-07-25, F.Osorio */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern void dupl_cols(int *, int *);

void
dupl_cols(int *order, int *cols)
{
  int n = *order, k = 0;

  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      k++;
      cols[i + j * n] = cols[j + i * n] = k;
    }
  }
}
