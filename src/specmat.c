/* $ID: specmat.c, last updated 2020-08-01, F.Osorio */

#include "base.h"
#include "specmat.h"

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
