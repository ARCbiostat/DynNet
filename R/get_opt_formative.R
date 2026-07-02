
#' Recover Lambda-scale parameters and formative weights after optimization
#'
#' @param temp list returned by optimizer (contains temp$b = optimized parameters)
#' @param paras_block_dim vector giving the sizes of parameter blocks
#' @param mapping vector: outcome → latent process (endogenous)
#' @param mapping2 vector: outcome → latent process (exogenous)
#'
#' @return temp list with parameters re-expressed in Lambda-space
#' @export

get_opt_formative <- function(b, paras_block_dim, mapping, mapping2, nRE, varcovRE.format ) {
  
  ## ---- (0) Extract Ω-space optimized parameters ----
  opt <- b
  pos <- c(1, cumsum(paras_block_dim))
  
  ## ---- (1) Mapping Ω → Λ ----
  mappingLP2LP1 <- pmin(table(mapping, mapping2), 1)
  nL <- ncol(mappingLP2LP1)
  nD <- nrow(mappingLP2LP1)
  
  mapping_vec <- apply(mappingLP2LP1, 2, function(x) which(x != 0))
  
  ## ---- (2) Extract Ω parameters ----
  alpha_mu0_trans <- opt[pos[1]:pos[2]]
  alpha_mu_trans  <- opt[(pos[2]+1):pos[3]]
  alpha_D         <- opt[(pos[3]+1):pos[4]]
  vec_alpha_ij    <- opt[(pos[4]+1):pos[5]]
  
  ## ---- (3) B_Omega ----
  if(varcovRE.format =="cholesky"){
    B_Omega <- DparChol(nRE,alpha_D)
  }
  if(varcovRE.format =="block"){
    B_Omega <- DparBlock(nD,(nb_RE-nD)/nD,alpha_D)
  }
  
  ## ---- (4) Baseline block (int0) ----
  B_u <- B_Omega[1:nL, 1:nL]
  
  ## ---- (5) RAW weights (somma = 1) ----
  weights_list <- vector("list", nD)
  
  for(d in 1:nD){
    idx_k <- which(mapping_vec==d)
    B_sub <- B_u[idx_k, idx_k, drop=FALSE]
    
    C_d <- rowSums(B_sub)
    w_raw <- solve(B_sub, C_d)
    w <- w_raw / sum(w_raw)
    
    weights_list[[d]] <- w
  }
  
  ## ---- (6) Build W_raw ----
  W_raw <- matrix(0, nrow=nD, ncol=nL)
  for(d in 1:nD){
    W_raw[d, which(mapping_vec==d)] <- weights_list[[d]]
  }
  
  ## ---- (7) Build We (from W_raw) ----
  W <- W_raw  # copy (verrà modificato)
  
  We <- matrix(0, nrow=sum(nrow(W)*(nRE/nL)), ncol=ncol(W)*(nRE/nL))
  We[1:nD,1:nL] <- W 
  
  for (d in 1:nD){
    for ( p in 1:((nRE/nL)-1)){
      col <- min(which(apply(We,2,sum)==0))
      row <- min(which(apply(We,1,sum)==0))
      w <- W[d,which(mapping_vec==d)]
      We[row,col:(col+length(w)-1)] <- w
    }
  }
  
  ## ---- (8) B_Lambda ----
  B_Lambda <- We %*% B_Omega %*% t(We)
  
  ## ===============================
  ## ✅ NORMALIZZAZIONE SOLO INT0
  ## ===============================
  
  scaling <- rep(1, nrow(B_Lambda))
  scaling[1:nD] <- 1 / sqrt(diag(B_Lambda)[1:nD])
  D <- diag(scaling)
  
  # covarianza normalizzata
  B_Lambda <- D %*% B_Lambda %*% D
  
  # pesi scalati (solo sulle righe int0)
  W <- D[1:nD, 1:nD] %*% W
  
  ## ---- (8b) rebuild We coerente ----
  We_new <- matrix(0, nrow=nrow(We), ncol=ncol(We))
  We_new[1:nD,1:nL] <- W 
  
  for (d in 1:nD){
    for ( p in 1:((nRE/nL)-1)){
      col <- min(which(apply(We_new,2,sum)==0))
      row <- min(which(apply(We_new,1,sum)==0))
      w <- W[d,which(mapping_vec==d)]
      We_new[row,col:(col+length(w)-1)] <- w
    }
  }
  We <- We_new
  
  ## CHECK
  # print(diag(B_Lambda)[1:nD])
  
  ## ===============================
  ## ✅ alpha con W SCALATI
  ## ===============================
  
  n_mu0 <- length(alpha_mu0_trans)/nL
  alpha_mu0 <- numeric(n_mu0 * nD)
  
  for(d in 1:nD){
    idx_k <- which(mapping_vec == d)
    w_d   <- W[d, idx_k]
    
    for(p in 1:n_mu0){
      k_idx <- idx_k + (p-1)*nL
      alpha_mu0[(d-1)*n_mu0 + p] <- sum(w_d * alpha_mu0_trans[k_idx])
    }
  }
  
  n_mu <- length(alpha_mu_trans)/nL
  alpha_mu <- numeric(n_mu * nD)
  
  for(d in 1:nD){
    idx_k <- which(mapping_vec == d)
    w_d   <- W[d, idx_k]
    
    for(p in 1:n_mu){
      k_idx <- idx_k + (p-1)*nL
      alpha_mu[(d-1)*n_mu + p] <- sum(w_d * alpha_mu_trans[k_idx])
    }
  }
  
  ## ---- (9) A_Lambda ----
  A_Omega <- matrix(vec_alpha_ij, nrow = nL, ncol = nL, byrow = TRUE)
  A_Lambda <- W %*% A_Omega %*% t(W)   
  vec_alpha_ij <- as.numeric(t(A_Lambda))
  
  ## ===============================
  ## ✅ SALVATAGGIO PESI
  ## ===============================
  
  para_weights_raw    <- as.numeric(t(W_raw))  # interpretabili
  para_weights_scaled <- as.numeric(t(W))      # usati nel modello
  
  ## ---- (10) altri parametri ----
  paraB <- ifelse(pos[5]==pos[6],numeric(0),opt[(pos[5]+1):pos[6]])
  paraSig <- opt[(pos[6]+1):pos[7]]
  ParaTransformY <- opt[(pos[8]+1):pos[9]]
  para_surv <- ifelse(pos[9]==pos[10],numeric(0),opt[(pos[9]+1):pos[10]])
  
  ## ---- (11) output ----
  b_inv <- c(
    alpha_mu0,
    alpha_mu,
    as.numeric( B_Lambda[ lower.tri(B_Lambda, diag=TRUE) ] ),
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv,
    para_weights_scaled   
  )
  
  ## restituisci anche raw per interpretazione
  return(list(
    b = b_inv,
    weights_raw = para_weights_raw,
    weights_scaled = para_weights_scaled
  ))
}
