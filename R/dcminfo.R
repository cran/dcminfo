#' @include dcminfo.R
#' @title The information matrix for diagnostic classification models
#' (or cognitive diagnostic models)
#'
#' @description This function is used to estimate the information matrix for
#' diagnostic classification models (DCMs; Rupp, Templin, & Henson, 2010) or cognitive diagnostic models,
#' such as the the observed information matrix,
#' the empirical cross-product information matrix and the sandwich-type
#' covariance matrix that can be used to estimate the asymptotic covariance
#' matrix or the model parameter standard errors.
#' @param dat A  \eqn{N \times J} binary data \code{matrix} consisting of the
#' responses of \eqn{N} examinees to \eqn{J} items.
#' @param delta A \code{list} of item parameter estimates.
#' @param attr_probs A \code{vector} of the estimated attribute mastery profile probability.
#' @param q_matrix A \eqn{J \times K}  \code{matrix} (Q-matrix) defines which attributes are measured by which items.
#' @param Mj A \code{list} of the design matrices and labels for each item (see de la Torre, 2011).
#' @param Aj A \code{list} of the possible combinations of the required attributes for each item (see de la Torre, 2011).
#' @param attr_mast_patt A \eqn{L \times K} binary \code{matrix} defines the attribute mastery profiles. see \code{amps}.
#' @param linkfct Type of the link function for the DCMs. It can be "logit", "identity", or "log".
#' In the current version, only the "logit" link function is now available.
#' @param info_type The returned information (or covariance) matrix type. The
#'  It can be "XPD" (the empirical cross-product information matrix),
#' "Obs"(the observed information matrix), or "Sw" (the sandwich-type covariance matrix).
#' The default is "Sw".
#'
#' @return A \code{matrix} giving information, or covariance matrix.
#'
#' @author {Yanlou Liu, Qufu Normal University, \email{liuyanlou@@163.com} \cr Tao Xin, Beijing Normal University}
#' @export
#' @examples
#'
#'#Example 1.
#'#The sandwich-type covariance matrix, the empirical cross-product information matrix,
#'#and the observed information matrix for the DINA model
#'
#'simresp <- sim_DINA_N1000$simresp
#'head(simresp)
#'
#'simdelta <- sim_DINA_N1000$simdelta
#'simdelta
#'
#'simqmatrix <- sim_DINA_N1000$simqmatrix
#'simqmatrix
#'
#'simAj <- sim_DINA_N1000$simAj
#'simAj
#'
#'simMj <- sim_DINA_N1000$simMj
#'simMj
#'
#'simAttrProbs <- sim_DINA_N1000$simAttrProbs
#'simAttrProbs
#'
#'# The number of the item parameters
#'N_delta <- length(unlist(simdelta))
#'N_delta
#'
#'
#' #Example 1.1 The sandwich-type covariance matrix
#'
#'Sw_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                  q_matrix=simqmatrix, Mj=simMj, Aj=simAj)
#'
#'
#'Sw_se_delta <- sqrt(diag(Sw_res))[1:N_delta]
#'
#'Sw_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=Sw_se_delta)
#'Sw_est_delta_se
#'
#'
#' #Example 1.2 The empirical cross-product information matrix
#'
#'XPD_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                   q_matrix=simqmatrix, Mj=simMj, Aj=simAj, info_type = "XPD")
#'
#'# Calculate the covariance matrix of the model parameters based on the XPD matrix
#'inv_XPD_res <- solve(XPD_res)
#'
#'XPD_se_delta <- sqrt(diag(inv_XPD_res))[1:N_delta]
#'XPD_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=XPD_se_delta)
#'XPD_est_delta_se
#'
#'
#' #Example 1.3 The observed information matrix
#'
#'Obs_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                   q_matrix=simqmatrix, Mj=simMj, Aj=simAj, info_type = "Obs")
#'
#'# Calculate the covariance matrix of the model parameters based on the Obs matrix
#'inv_Obs_res <- solve(Obs_res)
#'
#'Obs_se_delta <- sqrt(diag(inv_Obs_res))[1:N_delta]
#'Obs_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=Obs_se_delta)
#'Obs_est_delta_se
#'
#'
#'# Example 2.
#'#The sandwich-type covariance matrix, the empirical cross-product information matrix,
#'#and the observed information matrix for the C-RUM
#'
#'simresp <- sim_CRUM_N1000$simresp
#'head(simresp)
#'
#'simdelta <- sim_CRUM_N1000$simdelta
#'simdelta
#'
#'simqmatrix <- sim_CRUM_N1000$simqmatrix
#'simqmatrix
#'
#'simAj <- sim_CRUM_N1000$simAj
#'simAj
#'
#'simMj <- sim_CRUM_N1000$simMj
#'simMj
#'
#'simAttrProbs <- sim_CRUM_N1000$simAttrProbs
#'simAttrProbs
#'
#'# The number of the item parameters
#'N_delta <- length(unlist(simdelta))
#'N_delta
#'
#'
#' #Example 2.1 The sandwich-type covariance matrix
#'
#'Sw_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                  q_matrix=simqmatrix, Mj=simMj, Aj=simAj)
#'
#'
#'Sw_se_delta <- sqrt(diag(Sw_res))[1:N_delta]
#'
#'Sw_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=Sw_se_delta)
#'Sw_est_delta_se
#'
#'
#' #Example 2.2 The empirical cross-product information matrix
#'
#'XPD_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                   q_matrix=simqmatrix, Mj=simMj, Aj=simAj, info_type = "XPD")
#'
#'# Calculate the covariance matrix of the model parameters based on the XPD matrix
#'inv_XPD_res <- solve(XPD_res)
#'
#'XPD_se_delta <- sqrt(diag(inv_XPD_res))[1:N_delta]
#'XPD_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=XPD_se_delta)
#'XPD_est_delta_se
#'
#'
#' #Example 2.3 The observed information matrix
#'
#'Obs_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#'                   q_matrix=simqmatrix, Mj=simMj, Aj=simAj, info_type = "Obs")
#'
#'# Calculate the covariance matrix of the model parameters based on the Obs matrix
#'inv_Obs_res <- solve(Obs_res)
#'
#'Obs_se_delta <- sqrt(diag(inv_Obs_res))[1:N_delta]
#'Obs_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=Obs_se_delta)
#'Obs_est_delta_se
#'
#'
#'#Example 3. User-specified attribute mastery patterns
#'
#'attr_mast_patt <- amps(q_matrix=simqmatrix)
#'
#'Sw_res <- dcminfo(dat=simresp, delta=simdelta, attr_probs=simAttrProbs,
#' attr_mast_patt = attr_mast_patt, q_matrix=simqmatrix, Mj=simMj, Aj=simAj)
#'
#'Sw_se_delta <- sqrt(diag(Sw_res))[1:N_delta]
#'
#'Sw_est_delta_se <- data.frame(delta_est= unlist(simdelta), se_delta=Sw_se_delta)
#'Sw_est_delta_se
#'
#'
#'#Example 4. Using the gdina function from the CDM package
#'library("CDM")
#'d1 <- CDM::gdina(data = sim_DINA_N1000$simresp, q.matrix = sim_DINA_N1000$simqmatrix,
#' maxit= 1000, rule="DINA", linkfct = "logit", calc.se=FALSE)
#'delta <- d1$delta
#'N_delta <- length(unlist(delta))
#'attr_probs <- d1$control$attr.prob[,1]
#'attr_mast_patt <- amps(q_matrix=simqmatrix)
#'Mj <- d1$Mj
#'Aj <- d1$Aj
#'
#'
#'#Example 4.1 The sandwich-type covariance matrix
#'
#'Sw_res <- dcminfo(dat=sim_DINA_N1000$simresp, delta=delta, attr_probs=attr_probs,
#'                  q_matrix=sim_DINA_N1000$simqmatrix, Mj=Mj, Aj=Aj)
#'
#'Sw_se_delta <- sqrt(diag(Sw_res))[1:N_delta]
#'
#'Sw_est_delta_se <- data.frame(delta_est= unlist(delta), se_delta=Sw_se_delta)
#'Sw_est_delta_se
#'
#'
#'#Example 4.2 The empirical cross-product information matrix
#'
#'XPD_res <- dcminfo(dat=sim_DINA_N1000$simresp, delta=delta, attr_probs=attr_probs,
#'                   q_matrix=simqmatrix, Mj=Mj, Aj=Aj, info_type = "XPD")
#'
#'# Calculate the covariance matrix of the model parameters based on the XPD matrix
#'inv_XPD_res <- solve(XPD_res)
#'
#'XPD_se_delta <- sqrt(diag(inv_XPD_res))[1:N_delta]
#'XPD_est_delta_se <- data.frame(delta_est= unlist(delta), se_delta=XPD_se_delta)
#'XPD_est_delta_se
#'
#'
#'#Example 4.3 The observed information matrix
#'
#'Obs_res <- dcminfo(dat=sim_DINA_N1000$simresp, delta=delta, attr_probs=attr_probs,
#'                   q_matrix=simqmatrix, Mj=Mj, Aj=Aj, info_type = "Obs")
#'
#'# Calculate the covariance matrix of the model parameters based on the Obs matrix
#'inv_Obs_res <- solve(Obs_res)
#'
#'Obs_se_delta <- sqrt(diag(inv_Obs_res))[1:N_delta]
#'Obs_est_delta_se <- data.frame(delta_est= unlist(delta), se_delta=Obs_se_delta)
#'Obs_est_delta_se
#'
#'
#' @references
#'
#' Liu, Y., Tian, W., & Xin, T. (2016). An Application of M2 Statistic to Evaluate the Fit of Cognitive Diagnostic Models. \emph{Journal of Educational and Behavioral Statistics, 41}, 3-26.
#'
#' Liu, Y., Xin, T., Andersson, B. & Tian, W. (2017). Information Matrix Estimation Procedures for Cognitive Diagnostic Models.
#' \emph{under review}.
#'
#' de la Torre, J. (2011). The generalized DINA model framework. \emph{Psychometrika, 76}, 179-199.
#'
#'Rupp, A. A., Templin, J., & Henson, R. A. (2010). \emph{Diagnostic measurement: theory, methods, and applications}. New York, NY: Guilford.


dcminfo <- function(dat, delta, attr_probs, q_matrix, Mj, Aj, attr_mast_patt = NULL, linkfct = "logit", info_type = "Sw"){

  dat <- as.matrix(dat)
  q_matrix <- as.matrix(q_matrix)
  if (is.null(attr_mast_patt)){
    attr_mast_patt <- amps(q_matrix)
  }
  if (linkfct != "logit"){
    stop("Sorry, only the logit link function is now available.")
  }
  all_info_type <- c("XPD" , "Obs" , "Sw")
 if (!(info_type %in% all_info_type)){
   stop("Please specify a correct information matrix type, such as XPD, Obs, or Sw.")
 }

  infoXPD <- lcdmxpd(dat = dat , delta = delta,
	attr_probs = attr_probs, q_matrix = q_matrix, Mj = Mj, Aj = Aj,
	attr_mast_patt = attr_mast_patt, linkfct = linkfct)

	if (info_type == "XPD"){
	  return(infoXPD)
	}

	if(info_type == "Obs") {
	infop2 <- lcdmobsp2(dat = dat , delta = delta,
	attr_probs = attr_probs, q_matrix = q_matrix, Mj = Mj, Aj = Aj,
	attr_mast_patt = attr_mast_patt, linkfct = linkfct)
	info <- infoXPD - infop2
	}

	if(info_type == "Sw") {
	infop2 <- lcdmobsp2(dat = dat , delta = delta,
	attr_probs = attr_probs, q_matrix = q_matrix, Mj = Mj, Aj = Aj,
	attr_mast_patt = attr_mast_patt, linkfct = linkfct)
	infoobs <- infoXPD - infop2
	invobsinfo <- SolveMat(infoobs)
	info <- invobsinfo %*% infoXPD %*% invobsinfo
	}

	return(info)

}
