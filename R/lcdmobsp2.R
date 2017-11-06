
lcdmobsp2 <- function(dat, delta, attr_probs, q_matrix, Mj, Aj, attr_mast_patt, linkfct){
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

  Npars <- length(unlist(delta)) + L -1

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

  bmpO <- as.vector(n_obs_patts/N)
  bmpiO  <- rowSums(p_xi_al_pal)
  pdptemp <- bmpO / (bmpiO)
  pdp <-  outer(pdptemp,rep(1,L) )


Info_ipar <- matrix(NA,Npars-L+1,Npars-L+1)

for (jj in 1:J){ # jj <- 1
  ajj <- ( attr_ind[[jj]] )
  mjjj <- Mj[[jj]][[1]][ ajj , ]
  M1 <- ncol(mjjj)
  itmself <- outer(seq(1:M1),rep(1,2))
  apairs0 <- t(combn(M1,2 ))
  parpairs <- rbind(itmself,apairs0)
  nparpairs <- nrow(parpairs)

  for (kk in 1:nparpairs){# kk <- 2
    mappair <- t(mjjj)[parpairs[kk,1],] * t(mjjj)[parpairs[kk,2],]
    mapplus <- outer( rep(1,n_uniq_resp), mappair )
    x1 <- outer( obs_resp_patt[,jj] , rep(1,L) )
    pjjjM <- outer( rep(1,n_uniq_resp) , p_xj_al[jj,] ) + 10^(-20)
    pg1 <- mapplus* pdp * p_xi_al_pal * (( x1 - pjjjM)^2 - pjjjM * (1 - pjjjM))
    indI1 <- Mj_index[jj,2]+parpairs[kk,1]-1
    indI2 <- Mj_index[jj,2]+parpairs[kk,2]-1
    Info_ipar[indI2,indI1] <-Info_ipar[indI1,indI2] <- sum( pg1 )
    }
}

itmpairscmb  <- t(combn(J,2 ))
itmpchoose <- nrow(itmpairscmb)
for (itp in 1:itmpchoose){ # itp <- 1
  itm1 <- itmpairscmb[itp,1]
  itm2 <- itmpairscmb[itp,2]

  ajj1 <- ( attr_ind[[itm1]] )
  mjjj1 <- Mj[[itm1]][[1]][ ajj1 , ]
  ajj2 <- ( attr_ind[[itm2]] )
  mjjj2 <- Mj[[itm2]][[1]][ ajj2 , ]
  M11 <- ncol(mjjj1)
  M22 <- ncol(mjjj2)

  for(k1 in 1:M11){
    rowind <- Mj_index[itm1,2] + k1 -1
    for(k2 in 1:M22){
      colind <- Mj_index[itm2,2] + k2 -1
      mappair <- t(mjjj1)[k1,] * t(mjjj2)[k2,]
      mapplus <- outer( rep(1,n_uniq_resp), mappair )
      x1 <- outer( obs_resp_patt[,itm1] , rep(1,L) )
      pjjjM <- outer( rep(1,n_uniq_resp) , p_xj_al[itm1,] ) + 10^(-20)
      x2 <- outer( obs_resp_patt[,itm2] , rep(1,L) )
      pjjjM2 <- outer( rep(1,n_uniq_resp) , p_xj_al[itm2,] ) + 10^(-20)
      pg1 <- mapplus* pdp * p_xi_al_pal * ( x1 - pjjjM) * ( x2 - pjjjM2)
      Info_ipar[colind,rowind] <- Info_ipar[rowind,colind] <- sum( pg1 )
    }
  }
}




Info_alpha <- matrix(0,L-1,L-1)

for(ll in 1:(L-1)){# ll <- 1
  DBeta <- GetOrder2ParBeta_l(Prior=attr_probs,l1=ll,l2=ll)
  temp <- outer( rep(1,n_uniq_resp), DBeta)
  Info_alpha[ll,ll] <- sum(pdp *p_xi_al * temp)
}


lpairs <- t(combn((L-1),2))
nlpairs <- choose((L-1),2)

for (numl in 1:nlpairs){ # numl <- 2
  l1 <- lpairs[numl,1]
  l2 <- lpairs[numl,2]
  DBeta <- GetOrder2ParBeta_l(Prior=attr_probs,l1=l1,l2=l2)
  temp <- outer( rep(1,n_uniq_resp), DBeta)
  Info_alpha[l2,l1] <- Info_alpha[l1,l2] <- sum(pdp *p_xi_al * temp)
}


Infoiparpalpha <- matrix(0,Npars-L+1,L-1)
ind <- 1
for (jj in 1:J){ # jj <- 15
  ajj <- ( attr_ind[[jj]] )
  mjjj <- Mj[[jj]][[1]][ ajj , ]
  M1 <- ncol(mjjj)
  for (kk in 1:M1){# kk <- 1
    mappair <- t(mjjj)[kk,]
    mapplus <- outer( rep(1,n_uniq_resp), mappair )
    x1 <- outer( obs_resp_patt[,jj] , rep(1,L) )
    pjjjM <- outer( rep(1,n_uniq_resp) , p_xj_al[jj,] ) + 10^(-20)
    for(ll in 1:(L-1)){# ll <- 1
      DBeta <- GetParBeta_l(Prior=attr_probs,ll=ll)
      temp <- outer( rep(1,n_uniq_resp), DBeta)
      pg1 <- mapplus * pdp * temp * p_xi_al * ( x1 - pjjjM)
      Infoiparpalpha[ind,ll] <- sum(pg1)
    }
    ind <- ind + 1
  }
}

Infominus <- rbind(cbind(Info_ipar, Infoiparpalpha) , cbind(t(Infoiparpalpha),Info_alpha))

resInfo <- Infominus * N
resInfo
}
