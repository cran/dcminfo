#' @include dcminfo.R
#' @title  Generate the possible attribute mastery profile
#'
#' @description This function is used to generate the possible attribute mastery profile.
#' @param ind_zero_probs a  \code{vector} of integers indicating which attribute mastery profiles have zero probability.
#' The default is NULL (none of the attribute mastery profiles has zero probability).
#' @param q_matrix A \eqn{J \times K} \code{matrix} (Q-matrix) defines which attributes are measured by which items.
#' @return A \code{matrix} giving the generated attribute mastery profiles.
#'
#'@author {Yanlou Liu, Qufu Normal University, \email{liuyanlou@@163.com} \cr Tao Xin, Beijing Normal University}
#'
#' @export
#' @examples
#'
#'# Example 1.
#'simqmatrix <- sim_DINA_N1000$simqmatrix
#'simqmatrix
#' attr_mast_patt <- amps(q_matrix=simqmatrix)
#' attr_mast_patt
#'
#'# Example 2.
#' ind_zero_probs <- c(3,7)
#'attr_mast_patt <- amps(q_matrix=simqmatrix, ind_zero_probs=ind_zero_probs)
#'attr_mast_patt
#'
#'
#' @references
#'George, A. C., Robitzsch, A., Kiefer, T., Gross, J., & Uenlue, A. (2016). The R Package CDM for cognitive diagnosis models. \emph{Journal of Statistical Software}, 74(2), 1-24. doi:10.18637/jss.v074.i02
#'
#'Ma, W. & de la Torre, J. (2017). GDINA: The generalized DINA model framework. \emph{R package version 1.4.2}. Retrived from https://CRAN.R-project.org/package=GDINA
#'
#'Robitzsch, A., Kiefer, T., George, A. C., & Uenlue, A. (2017). \emph{CDM: Cognitive diagnosis modeling. R package version 5.9-27}. Retrived from https://CRAN.R-project.org/package=CDM
#'

amps <- function(q_matrix,ind_zero_probs=NULL){
  q_matrix <- as.matrix(q_matrix)
  K <- ncol(q_matrix)
  if (max(q_matrix)!=1){
    stop("Only binary attributes available.")
  }
  if (K <= 1){
    stop("The number of attributes must be larger than one.")
  }
  n_k_list <- list()
  for (kk in 1:K){ n_k_list[[kk]] <- c(0,1) }

  amps_all <- as.matrix( expand.grid( n_k_list ) )
  if(!is.null(ind_zero_probs)){
    res_amp <- amps_all[-ind_zero_probs,]
  } else{
    res_amp <- amps_all
  }
  return(res_amp)
}
