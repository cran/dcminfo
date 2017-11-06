#' Simulated observed response data, Q-matrix and model parameters (the deterministic inputs, noisy "and" gate model)
#'
#' Artificial Q-matrix for 14 items with 3 attributes.
#' A list of observed response data for 14 items,
#' that was generated using the the deterministic inputs, noisy "and" gate model (DINA; Haertel, 1989; Junker & Sijtsma, 2001) .
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
#' Haertel, E. H. (1989). Using restricted latent class models to map the skill structure of achievement items.
#' \emph{Journal of Educational Measurement, 26}, 301-321.
#'
#' Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with few assumptions, and connections with nonparametric
#' item response theory. \emph{Applied Psychological Measurement, 25}, 258-272.
#'
"sim_DINA_N1000"
