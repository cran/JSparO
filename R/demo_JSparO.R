#' @title demo_JSparO - The demo of JSparO package
#' @description This is the main function of JSparO aimed to solve the low-order regularization models with \eqn{l_{p,q}} norm.
#' @param A Gene expression data of transcriptome factors (i.e. feature matrix in machine learning).
#' The dimension of A is m * n.
#' @param B Gene expression data of target genes (i.e. observation matrix in machine learning).
#' The dimension of B is m * t.
#' @param X Gene expression data of Chromatin immunoprecipitation or other matrix
#' (i.e. initial iterative point in machine learning). The dimension of X is n * t.
#' @param s joint sparsity level
#' @param p value for \eqn{l_{p,q}} norm (i.e., p = 1 or 2)
#' @param q value for \eqn{l_{p,q}} norm (i.e., 0 <= q <= 1)
#' @param maxIter maximum iteration
#' @return The solution of proximal gradient method with \eqn{l_{p,q}} regularizer.
#' @author Xinlin Hu <thompson-xinlin.hu@connect.polyu.hk>
#'
#' Yaohua Hu <mayhhu@szu.edu.cn>
#' @details The demo_JSparO function is used to solve joint sparse optimization problem via different algorithms.
#' Based on \eqn{l_{p,q}} norm, functions with different p and q are implemented to solve the problem:
#' \deqn{\min \|AX-B\|_F^2 + \lambda \|X\|_{p,q}}
#' to obtain s-joint sparse solution.
#' @export demo_JSparO
#' @examples
#' m <- 256; n <- 1024; t <- 5; maxIter0 <- 50
#' A0 <- matrix(rnorm(m * n), nrow = m, ncol = n)
#' B0 <- matrix(rnorm(m * t), nrow = m, ncol = t)
#' X0 <- matrix(0, nrow = n, ncol = t)
#' res_JSparO <- demo_JSparO(A0, B0, X0, s = 10, p = 2, q = 'half', maxIter = maxIter0)
#'
demo_JSparO <- function(A, B, X, s, p, q, maxIter = 200){
  # Normalization
  NoA <- norm(A, '2'); A <- A/NoA; B <- B/NoA

  if(q == '0'){
    alg_name <- paste0('L', p, 'HardThr(A, B, X, s, maxIter)')
  }else if(q == 'half'){
    alg_name <- paste0('L', p, 'HalfThr(A, B, X, s, maxIter)')
  }else if(q == 'twothirds'){ #
    alg_name <- paste0('L', p, 'twothirdThr(A, B, X, s, maxIter)')
  }else if(q == '1'){ # L2Soft
    alg_name <- paste0('L', p, 'SoftThr(A, B, X, s, maxIter)')
  }else{ # Newton method
    alg_name <- paste0('L', p, 'NewtonThr(A, B, X, s, q, maxIter)')
  }

  res <- eval(parse(text = alg_name))
  message('Success.\n')
  return(res)
}
