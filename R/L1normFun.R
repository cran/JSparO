#' @title L1normFun
#' @description The function aims to compute the \eqn{l_1} norm.
#' @param x vector
#' @return The \eqn{l_1} norm of vector x
#' @author Xinlin Hu <thompson-xinlin.hu@connect.polyu.hk>
#'
#' Yaohua Hu <mayhhu@szu.edu.cn>
#' @details The L1normFun aims to compute the \eqn{l_1} norm: \eqn{\sum_i^n |x_i|}
#'
L1normFun <- function(x){
  return(sum(abs(x)))
}
