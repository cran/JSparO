#' @title L2twothirdsThr - Iterative Thresholding Algorithm based on \eqn{l_{2,2/3}} norm
#' @description The function aims to solve \eqn{l_{2,2/3}} regularized least squares.
#' @param A Gene expression data of transcriptome factors (i.e. feature matrix in machine learning).
#' The dimension of A is m * n.
#' @param B Gene expression data of target genes (i.e. observation matrix in machine learning).
#' The dimension of B is m * t.
#' @param X Gene expression data of Chromatin immunoprecipitation or other matrix
#' (i.e. initial iterative point in machine learning). The dimension of X is n * t.
#' @param s joint sparsity level
#' @param maxIter maximum iteration
#' @return The solution of proximal gradient method with \eqn{l_{2,2/3}} regularizer.
#' @author Xinlin Hu <thompson-xinlin.hu@connect.polyu.hk>
#'
#' Yaohua Hu <mayhhu@szu.edu.cn>
#' @details The L2twothirdsThr function aims to solve the problem:
#' \deqn{\min \|AX-B\|_F^2 + \lambda \|X\|_{2,2/3}}
#' to obtain s-joint sparse solution.
#' @export L2twothirdsThr
#' @examples
#' m <- 256; n <- 1024; t <- 5; maxIter0 <- 50
#' A0 <- matrix(rnorm(m * n), nrow = m, ncol = n)
#' B0 <- matrix(rnorm(m * t), nrow = m, ncol = t)
#' X0 <- matrix(0, nrow = n, ncol = t)
#' NoA <- norm(A0, '2'); A0 <- A0/NoA; B0 <- B0/NoA
#' res_L2twothirds <- L2twothirdsThr(A0, B0, X0, s = 10, maxIter = maxIter0)
#'
L2twothirdsThr <- function(A, B, X, s, maxIter = 200){
  # Initialization
  v <- 0.5 # stepsize
  Va1 <- (0.5)^(4/3) * 1.5 / v
  Bu1 <- 2 * v * t(A) %*% B
  Bu2 <- 2 * v * t(A) %*% A

  for(k in 1:maxIter){
    # Gradient descent
    Bu <- X + Bu1 - Bu2 %*% X

    # L2-2/3 threhsolding operator
    normBu <- apply(Bu, 1, norm, '2') # L2 norm of each row
    criterion <- sort(normBu, decreasing = T)[s+1]
    lambda <- Va1 * criterion^(4/3)
    q <- 2 * lambda * v

    # Consider what if s-th largest group is not the only one
    if(criterion == sort(normBu, decreasing = T)[s]){
      ind <- which(normBu >= criterion)
    }else{
      ind <- which(normBu > criterion)
    }

    Xnew <- matrix(0, nrow = nrow(Bu), ncol = ncol(Bu))
    for(i in ind){
      rowDa <- Bu[i, ]
      normRow <- norm(rowDa, '2')
      phi <- acosh(27 * normRow^2 / (16 * q^(1.5)))
      aa <- 2 * q^0.25 * (cosh(phi/3))^(0.5) / (sqrt(3))
      eta <- 3 * (aa^(1.5) + sqrt(2 * normRow - aa^3))^4

      # Update matrix
      Xnew[i, ] <- (eta/(32 * lambda * v * aa^2 + eta)) * rowDa
    }

    if(norm(Xnew - X, 'f') <= 1e-6){ break }
    # Update and report
    X <- Xnew
    # message(paste('The', k, '-th iteration completes.'))
  }

  return(Xnew)
}
