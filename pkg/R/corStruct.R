## ID: corStruct.R, last updated 2022-07-01, F.Osorio

corAR1 <- function(rho, p = 2)
{ ## form an AR(1) correlation matrix
  if ((rho > 1) || (rho < -1))
    stop("'rho' must be in the interval (-1,1)")
  if (p < 2)
    stop("dimension 'p' must be greater or equal than 2")

  Cor <- matrix(0, nrow = p, ncol = p)
  storage.mode(Cor) <- "double"
  o <- .C("cor_AR1",
          Cor = Cor,
          p = as.integer(p),
          rho = as.double(rho))$Cor
  o
}

corCS <- function(rho, p = 2)
{ ## form an equicorrelation (compound symmetry) matrix
  if ((rho > 1) || (rho < -1))
    stop("'rho' must be in the interval (-1,1)")
  if (p < 2)
    stop("dimension 'p' must be greater or equal than 2")

  Cor <- matrix(0, nrow = p, ncol = p)
  storage.mode(Cor) <- "double"
  o <- .C("cor_CS",
          Cor = Cor,
          p = as.integer(p),
          rho = as.double(rho))$Cor
  o
}
