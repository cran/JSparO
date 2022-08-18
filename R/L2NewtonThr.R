#' @title L2NewtonThr - Iterative Thresholding Algorithm based on \eqn{l_{2,q}} norm with Newton method
#' @description The function aims to solve \eqn{l_{2,q}} regularized least squares, where the proximal optimization subproblems will be solved by Newton method.
#' @param A Gene expression data of transcriptome factors (i.e. feature matrix in machine learning).
#' The dimension of A is m * n.
#' @param B Gene expression data of target genes (i.e. observation matrix in machine learning).
#' The dimension of B is m * t.
#' @param X Gene expression data of Chromatin immunoprecipitation or other matrix
#' (i.e. initial iterative point in machine learning). The dimension of X is n * t.
#' @param s joint sparsity level
#' @param q value for \eqn{l_{2,q}} norm (i.e., 0 < q < 1)
#' @param maxIter maximum iteration
#' @param innMaxIter maximum iteration in Newton step
#' @param innEps criterion to stop inner iteration
#' @return The solution of proximal gradient method with \eqn{l_{2,q}} regularizer.
#' @author Xinlin Hu <thompson-xinlin.hu@connect.polyu.hk>
#'
#' Yaohua Hu <mayhhu@szu.edu.cn>
#' @details The L2NewtonThr function aims to solve the problem:
#' \deqn{\min \|AX-B\|_F^2 + \lambda \|X\|_{2,q}}
#' to obtain s-joint sparse solution.
#' @import pracma
#' @export L2NewtonThr
#' @examples
#' m <- 256; n <- 1024; t <- 5; maxIter0 <- 50
#' A0 <- matrix(rnorm(m * n), nrow = m, ncol = n)
#' B0 <- matrix(rnorm(m * t), nrow = m, ncol = t)
#' X0 <- matrix(0, nrow = n, ncol = t)
#' NoA <- norm(A0, '2'); A0 <- A0/NoA; B0 <- B0/NoA
#' res_L2q <- L2NewtonThr(A0, B0, X0, s = 10, q = 0.2, maxIter = maxIter0)
#'
L2NewtonThr <- function(A, B, X, s, q, maxIter = 200, innMaxIter = 30, innEps = 1e-6){
  # Initialization
  v <- 0.5 # stepsize
  Bu1 <- 2 * v * t(A) %*% B
  Bu2 <- 2 * v * t(A) %*% A
  n <- ncol(B)
  I <- diag(1, n, n)

  for(k in 1:maxIter){
    # Gradient descent
    Bu <- X + Bu1 - Bu2 %*% X

    # Newton method to solve subproblem
    normBu <- apply(Bu, 1, norm, '2') # L2 norm of each row
    ind <- order(normBu, decreasing = T)[1:s]
    criterion <- sort(normBu, decreasing = T)[s+1]
    lambda <- criterion^(2-q) / v

    Xnew <- matrix(0, nrow = nrow(Bu), ncol = ncol(Bu))
    for(i in ind){
      rowDa <- Bu[i, ]
      rowDaTemp <- rowDa

      # Newton method
      for(j in 1:innMaxIter){
        normTemp <- norm(rowDaTemp, '2')
        H <- lambda * q * rowDaTemp + (rowDaTemp - rowDa) * normTemp^(2 - q) / v
        DH <- lambda * q * I + normTemp^(-q) * (normTemp^2 * I + (2 - q) * (rowDaTemp - rowDa) %*% t(rowDaTemp)) / v

        # rowDaTemp <- rowDaTemp - solve(DH) %*% H
        rowDaTemp <- tryCatch({
          rowDaTemp - solve(DH) %*% H
        }, error = function(cond){
          rowDaTemp - pinv(DH) %*% H
        })
        HTemp <- lambda * q * rowDaTemp + (rowDaTemp - rowDa) * norm(rowDaTemp, '2')^(2 - q) / v
        if(norm(HTemp, '1') < innEps){ break }
      }

      # Update matrix
      Xnew[i, ] <- rowDaTemp
    }

    if(norm(Xnew - X, 'f') <= 1e-6){ break }
    # Update and report
    X <- Xnew
    # message(paste('The', k, '-th iteration completes.'))
  }

  return(Xnew)
}
