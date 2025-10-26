## ID: floyd.R, last updated 2025-10-15, F.Osorio

floyd <- function(x)
{ ## find the shortest paths in a directed graph
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")

  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  if (n != p)
    stop("argument x is not a square matrix")

  BIG <- .Machine$double.xmax
  x[is.na(x)] <- BIG
  x[is.infinite(x)] <- BIG
  storage.mode(x) <- "double"

  paths <- matrix(0, nrow = n, ncol = n)
  storage.mode(paths) <- "integer"
  
  # initializing
  z <- .Fortran("floyd_init",
                costs = x, 
                ldc = as.integer(n), 
                paths = paths, 
                ldp = as.integer(n), 
                n = as.integer(n))[c("costs","paths")]
  
  # Floyd-Warshall iterations
  z <- .Fortran("floyd_warshall",
                costs = z$costs, 
                ldc = as.integer(n), 
                paths = z$paths, 
                ldp = as.integer(n), 
                n = as.integer(n))[c("costs","paths")]

  # output
  which <- z$costs == BIG
  if (any(which)) {
    z$costs[which] <- Inf
    z$paths[which] <- 0
  }
  z <- list(costs = z$costs, shortest.paths = z$paths)
  z
}
