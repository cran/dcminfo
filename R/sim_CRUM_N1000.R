#' Simulated observed response data, Q-matrix and model parameters (compensatory reparameterized unified model)
#'
#' Artificial Q-matrix for 14 items with 3 attributes.
#' A list of observed response data for 14 items,
#' that was generated using the compensatory reparameterized unified model (C-RUM; Hartz, 2002).
#'
#' @format A list of observed response data Q-matrix, model parameters for 14 items with components:
#' \describe{
#' \item{\code{simresp}}{simulated responses data matrix of 1000 examinees response to 14 items.}
#' \item{\code{simdelta}}{A list of simulated item parameters for 14 items.}
#' \item{\code{simqmatrix}}{Artificial Q-matrix specifies the relationship between 14 items and 3 attributes.}
#' \item{\code{simAj}}{A \code{list} of the possible combinations of the required attributes for 14 items (see de la Torre, 2011).}
#' \item{\code{simMj}}{A \code{list} of the design matrices and labels for 14 items (see de la Torre, 2011).}
#' \item{\code{simAttrProbs}}{A Simulated \code{vector} of the probabilities of attribute mastery profile of \eqn{2^3}.}
#' }
#' @references
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#' Hartz, S. M. (2002). A bayesian framework for the unified model for assessing cognitive abilities:
#' Blending theory with practicality (Unpublished doctoral dissertation). University of Illinois at Urbana-Champaign.
#'
"sim_CRUM_N1000"
