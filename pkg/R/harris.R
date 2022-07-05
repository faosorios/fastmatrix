## ID: harris.R, last updated 2022-07-01, F.Osorio

harris.test <-
function(x, test = "Wald")
{
  z <- as.matrix(x)
  if (!is.numeric(z))
    stop("harris applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]

  ## sample covariance matrix under normality
  S <- cov(z)

  ## contrast matrix
  H <- cbind(-1, diag(p - 1))

  switch(test,
    "Wald" = {
      s <- diag(S)
      S2 <- hadamard(S)
      G <- tcrossprod(H %*% S2, H)
      h <- H %*% s
      g <- solve(G, h)
      stat <- .5 * n * sum(g * h)
      names(stat) <- "Wald"
      method <- "Wald test"
    },
    "log" = {
      v <- log(diag(S))
      R <- cov2cor(S)
      R2 <- hadamard(R)
      G <- tcrossprod(H %*% R2, H)
      h <- H %*% v
      g <- solve(G, h)
      stat <- .5 * n * sum(g * h)
      names(stat) <- "Wald"
      method <- "log-transformation Wald test"
    },
    "robust" = {
      s <- diag(S)
      xbar <- apply(z, 2, mean)
      Psi <- matrix(0, nrow = p, ncol = p)
      storage.mode(Psi) <- "double"
      storage.mode(z) <- "double"
      Psi <- .C("cov4th",
              z = z,
              n = as.integer(n),
              p = as.integer(p),
              center = as.double(xbar),
              Psi = Psi)$Psi
      G <- tcrossprod(H %*% Psi, H)
      h <- H %*% s
      g <- solve(G, h)
      prod <- sum(g * h)
      stat <- n * prod / (1 - prod)
      names(stat) <- "Wald"
      method <- "robust Wald test"
    },
    "log-robust" = {
      s <- diag(S)
      v <- log(s)
      xbar <- apply(z, 2, mean)
      Q <- matrix(0, nrow = p, ncol = p)
      storage.mode(Q) <- "double"
      storage.mode(z) <- "double"
      Q <- .C("cov4th",
              z = z,
              n = as.integer(n),
              p = as.integer(p),
              center = as.double(xbar),
              Q = Q)$Q
      storage.mode(Q) <- "double"
      Q <- .C("Psi2Q",
              Q = Q,
              s = as.double(s),
              p = as.integer(p))$Q
      G <- tcrossprod(H %*% Q, H)
      h <- H %*% v
      g <- solve(G, h)
      stat <- n * sum(g * h)
      names(stat) <- "Wald"
      method <- "log-robust Wald test"
    },
    stop(paste("unimplemented option:", test))
  )
  df <- p - 1
  pval <- 1 - pchisq(stat, df = df)

  ## output object
  z <- list(statistic = stat, parameter = df, p.value = pval, estimate = S, method = method)
  class(z) <- "Harris.test"
  z
}

print.Harris.test <- function(x, digits = 4, ...)
{
  ## local function
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }

  cat("\n")
  cat(paste(x$method, "for equality of variances", sep = " "), "\n")
  cat("\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: true variances are not equal.\n")
  cat("\n")
  cat("sample estimate:\n")
  print.symmetric(x$estimate, digits = digits)
  invisible(x)
}
