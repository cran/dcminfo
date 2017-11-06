

GetOrder2ParBeta_l <- function(Prior,l1,l2){
  L <- length(Prior)
  D_value <- 1/Prior[L]
  expBeta <- D_value * Prior
  if (l1==l2){
    res <- -(expBeta * expBeta[l1])*(D_value - 2*expBeta[l1])/(D_value^3)
    res[l1] <- expBeta[l1]*(D_value - expBeta[l1])*(D_value - 2*expBeta[l1])/(D_value^3)
  } else {
    res <- 2*(expBeta * expBeta[l1]*expBeta[l2])/(D_value^3)
    res[l1] <- expBeta[l1]*expBeta[l2]*(-D_value + 2*expBeta[l1])/(D_value^3)
    res[l2] <- expBeta[l1]*expBeta[l2]*(-D_value + 2*expBeta[l2])/(D_value^3)
  }
  res
}



GetParBeta_l <- function(Prior,ll){
  L <- length(Prior)
  D_value <- 1/Prior[L]
  expBeta <- D_value * Prior
  expBeta_l <- expBeta[ll]
  res <- -(expBeta_l * expBeta)/(D_value^2)
  res[ll] <- expBeta_l*(D_value - expBeta_l)/(D_value^2)
  res
}
#' @importFrom methods is

SolveMat <- function(Info){
	eps2 <- 10^(-10)
	varmatTemp <- try( solve(Info + diag( eps2 , ncol(Info) )) )
	if ( is(varmatTemp , "try-error") ){
	cat( "Singular model parameter covariance matrix\n")
	varmatTemp <- NA}
	varmatTemp
}
