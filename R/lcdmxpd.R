#' @importFrom stats plogis
#' @importFrom utils combn
#'
#'
#'
lcdmxpd <- function(dat, delta, attr_probs, q_matrix, Mj, Aj, attr_mast_patt, linkfct){

datmat <- as.matrix(dat)
N <- nrow(datmat)
J <- ncol(datmat)
L <- length(attr_probs)

  n_obs_patts <- table(apply(datmat, 1, paste, collapse=""))
  obs_patt_list <- strsplit(rownames(n_obs_patts),split = NULL)
  n_uniq_resp <- length(n_obs_patts)
  obs_resp_patt <- matrix(NA,n_uniq_resp,J)
  for (orp in 1:n_uniq_resp){
    obs_resp_patt[orp,] <- as.numeric(obs_patt_list[[orp]])
  }

  Mj_index <- matrix( 0 , J , 3 )
  attr_ind <- list()
  for (jj in 1:J){
    attr_loc_j <- which( q_matrix[jj,] > 0 )
    amp_qmatrix <- 1 * (attr_mast_patt[ ,attr_loc_j] == q_matrix[ rep(jj,L) ,attr_loc_j])
    attr_ind_j <- rep(0,L)
    for (dd in 1:nrow(Aj[[jj]]) ){
      attr_ind_j[ which( rowMeans( amp_qmatrix == outer( rep(1,nrow(attr_mast_patt)) , Aj[[jj]][dd,] )  ) == 1) ] <- dd
    }
    attr_ind[[jj]] <- attr_ind_j
    Mj_index[jj,1] <- ncol( Mj[[jj]][[1]] )	#NOTE THIS
  }
  Mj_index[,3] <- cumsum( Mj_index[,1] )
  Mj_index[,2] <- c(1,Mj_index[-J,3] + 1 )

  p_xj_al <- matrix(0, nrow = J, ncol = L)
  for (jj in 1:J){
    ajj <- ( attr_ind[[jj]] )
    mjjj <- Mj[[jj]][[1]][ ajj , ]
    djj <- matrix( delta[[jj]] , L , length(delta[[jj]]) , byrow=TRUE )
    p_xj_al[jj,] <- rowSums( mjjj * djj )
    if (linkfct == "logit"){p_xj_al[jj,] <- plogis( p_xj_al[jj,] )}
  }

  p_xj_al[ p_xj_al < 0 ] <- 10^(-10)
  p_xj_al[ p_xj_al > 1] <- 1 - 10^(-10)
  p_xj0_al <- 1- p_xj_al

  p_xi_al <- matrix( 1, n_uniq_resp, L)

  for (lll in 1:L) {
    for(iii in 1:n_uniq_resp) {
      for(jjj in 1:J) {
        if(obs_resp_patt[iii,jjj]==1) { p_xi_al[iii,lll] <- p_xj_al[jjj, lll] * p_xi_al[iii,lll] }
        if(obs_resp_patt[iii,jjj]==0) { p_xi_al[iii,lll] <- p_xj0_al[jjj, lll]  * p_xi_al[iii,lll] }
      }
    }
  }

  p_xi_al_pal <- outer( rep(1,n_uniq_resp), attr_probs ) * p_xi_al

  Jac_ipar <- matrix(0,n_uniq_resp, length(unlist(delta)))

  for (jj in 1:J){ # jj <- 2
  ajj <- ( attr_ind[[jj]] )
  mjjj <- Mj[[jj]][[1]][ ajj , ]
  M1 <- ncol(mjjj)
  p_per_deltaj <- matrix( 0 , nrow=n_uniq_resp , ncol = M1 )
  x1 <- outer( obs_resp_patt[,jj] , rep(1,L) )
  pjjjM <- outer( rep(1,n_uniq_resp) , p_xj_al[jj,] ) + 10^(-20)
    for (kk in 1:M1){#  kk <- 1
      mapplus <- outer( rep(1,n_uniq_resp), t(mjjj)[kk,] )
	if(linkfct=="logit"){
		p_per_deltaj_jj <- mapplus * p_xi_al_pal * ( x1 - pjjjM) #/ ( pjjjM * ( 1 - pjjjM ) )
	}
      p_per_deltaj[,kk] <- rowSums( p_per_deltaj_jj )
    }
    Jac_ipar[,Mj_index[jj,2]:Mj_index[jj,3]] <- p_per_deltaj
  }
  Jac_pal <- matrix(0,n_uniq_resp,L-1)
  for(ll in 1:(L-1)){# ll <- 1
	DBeta <- GetParBeta_l(attr_probs,ll)
	temp <- outer( rep(1,n_uniq_resp), DBeta)
	Jac_pal[,ll] <- rowSums(p_xi_al * temp)
  }

  Jac <- cbind(Jac_ipar,Jac_pal)
  bmpO <- as.vector(n_obs_patts/N)
  bmpiO  <- rowSums(p_xi_al_pal)
  temp <- bmpO / (bmpiO^2)
  InfoO <- t(Jac) %*% diag(temp) %*% Jac
  return(InfoO*N);
}
