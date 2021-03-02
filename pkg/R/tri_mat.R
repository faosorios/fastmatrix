## ID: tri_mat.R, last updated 2020-10-19, F.Osorio

is.lower.tri <- function(x, diag = FALSE)
{ ## test if a matrix is lower triangular
  if (diag)
    ok <- all(x[row(x) >= col(x)] != 0)
  else
    ok <- all(x[row(x) > col(x)] != 0)
  ok
}

is.upper.tri <- function(x, diag = FALSE)
{ ## test if a matrix is upper triangular
  if (diag)
    ok <- all(x[row(x) <= col(x)] != 0)
  else
    ok <- all(x[row(x) < col(x)] != 0)
  ok
}
