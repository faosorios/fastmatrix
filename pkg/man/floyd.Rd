\name{floyd}
\alias{floyd}
\title{Find the shortest paths in a directed graph}
\usage{
floyd(x)
}
\arguments{
  \item{x}{ adjacency matrix of a directed graph. \eqn{(i,j)}-entry represents the 
    weight (cost) of the edge from \eqn{i} to \eqn{j} if one exists and \eqn{\infty}
    otherwise.}
}
\description{
  The Floyd-Warshall algorithm finds all shortest paths (if exist) in a directed graph.
}
\value{
  Returns a list with final costs and the shortest path between two nodes.
}
\references{
  Floyd, R.W. (1962). 
  Algorithm 97: Shortest Path. 
  \emph{Communications of the ACM} \bold{5} (6), 345.
}
\examples{
x <- matrix(c(0,3,Inf,5,2,0,Inf,4,Inf,1,0,Inf,Inf,Inf,2,0), nrow = 4, 
            ncol = 4, byrow = TRUE)
z <- floyd(x)
z
}
\keyword{array}
