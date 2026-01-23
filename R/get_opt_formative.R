
#' Recover Lambda-scale parameters and formative weights after optimization
#'
#' @param temp list returned by optimizer (contains temp$b = optimized parameters)
#' @param paras_block_dim vector giving the sizes of parameter blocks
#' @param mapping vector: outcome → latent process (endogenous)
#' @param mapping2 vector: outcome → latent process (exogenous)
#'
#' @return temp list with parameters re-expressed in Lambda-space
#' @export
get_opt_formative <- function(b, paras_block_dim, mapping, mapping2,nRE) {
  
  ## ---- (0) Extract Ω-space optimized parameters ----
  opt <- b
  pos <- c(1, cumsum(paras_block_dim))
  
  ## ------------ (1) Identify mapping Ω_k → Λ_d ----------------------
  # mappingLP2LP1 is a binary matrix: row = Λ_d, col = Ω_k
  mappingLP2LP1 <- pmin(table(mapping, mapping2), 1)
  nL <- ncol(mappingLP2LP1)  # number of exogenous processes
  nD <- nrow(mappingLP2LP1)  # number of endogenous processes
  
  # vector of Ω indices corresponding to each Λ
  mapping_vec <- apply(mappingLP2LP1, 2, function(x) which(x != 0))
  
  ## ---- (2) Extract Ω-space fixed effects (βk and γk) -----------
  alpha_mu0_trans <- opt[pos[1]:pos[2]]      # intercept terms for all Ω
  alpha_mu_trans  <- opt[(pos[2]+1):pos[3]]  # slope terms for all Ω
  alpha_D         <- opt[(pos[3]+1):pos[4]]  # Ω random effect covariance
  vec_alpha_ij    <- opt[(pos[4]+1):pos[5]]  # temporal effects (Ω-scale)
  
  ## ---- (3) Rebuild Ω-space random effect covariance matrix BΩ ----
  B_Omega <- DparChol(nRE,alpha_D)  # YOU ALREADY HAVE THIS FUNCTION
  
  ## ---- (4) Extract baseline Ω random effect covariance block ----
  # First nL elements correspond to u_k (random intercepts)
  # Make sure nL is known in this scope.
  B_u <- B_Omega[1:nL, 1:nL]
  
  ## ---- (5) Compute Λ_d(0) → needed for weight recovery ----
  # Λ_d(0) = X β_d + u_d
  # BUT β_d and u_d are not estimated directly.
  # They must be reconstructed from Ω-space parameters.
  
  
  # Compute Ω(0) = X0 * alpha_mu0_trans + u_k
  # But for the covariance-based weight formula, we use ONLY random effects:
  # Cov(Λ(0), Ω_k(0)) = Cov(u_d, u_k)
  
  ## ---- (6) Recover formative weights ------------------------------
  weights_list <- vector("list", nD)
  
  for(d in 1:nD){
    
    # which Ω processes contribute to Λ_d?
    idx_k <- which(mapping_vec==d)
    
    # extract the sub-covariance matrix among contributing Ω's
    B_sub <- B_u[idx_k, idx_k, drop=FALSE]
    
    # Cov(Λ_d(0), Ω_k(0)) = sum_l ω_l Cov(u_l, u_k)
    # But we don't know ω yet. However:
    # the regression coefficient formula gives:
    #   ω_k ∝ Cov(Λ_d(0), Ω_k(0))
    #
    # And Cov(Λ,Ω_k) is the k-th row-sum of B_sub * w, but Λ_d variance satisfies:
    #   Var(Λ_d) = w^T B_sub w
    
    # Simplest: set up the linear system
    # C_d = B_sub * w_d   with sum(w_d)=1
    #
    # But C_d is E[u_d * u_k] = Cov(u_d, u_k)
    # And Cov(u_d, u_k) = row-sum( B_sub )[k ]
    #
    # So C_d = rowSums(B_sub)
    
    C_d <- rowSums(B_sub)
    
    # Solve B_sub * w = C_d
    w_raw <- solve(B_sub, C_d)
    
    # Normalize to sum 1
    w <- w_raw / sum(w_raw)
    
    weights_list[[d]] <- w
  }
  
  ## ---- (7) Reconstruct Λ-scale fixed effects ----------------------
  # α_mu0(d) = sum_k ω_k α_mu0(k)  for each Λ_d
  
  n_mu0 <-  length(alpha_mu0_trans)/nL
  map_mu0 <-  rep(c(1:n_mu0),times=nL)
  mappingLP2LP1_vec_mu0 <-  rep(mapping_vec,each=n_mu0) 
  
  alpha_mu0 <- unlist(
    lapply(1:nD, function(d){lapply(1:n_mu0,function(p){
      idx_k <- which(mappingLP2LP1_vec_mu0==d & map_mu0==p)
      sum( weights_list[[d]] * alpha_mu0_trans[idx_k] )
    })})
  )
  
  print(alpha_mu0)
  
  n_mu <-  length(alpha_mu_trans)/nL
  map_mu <-  rep(c(1:n_mu),times=nL)
  mappingLP2LP1_vec_mu <-  rep(mapping_vec,each=n_mu) 
  
  alpha_mu <- unlist(
    lapply(1:nD, function(d){lapply(1:n_mu,function(p){
      idx_k <- which(mappingLP2LP1_vec_mu==d & map_mu==p)
      sum( weights_list[[d]] * alpha_mu_trans[idx_k] )
    })})
  )
  
  ## ---- (8) Reconstruct Λ-scale covariance -------------------------
  # B_Λ = W B_Ω Wᵀ
  # Build full W matrix
  W <- matrix(0, nrow=nD, ncol=nL)
  for(d in 1:nD){
    W[d, which(mapping_vec==d)] <- weights_list[[d]]
  }
  
  
  We <- matrix(0,nrow=sum(nrow(W)*(nRE/nL)),ncol=ncol(W)*(nRE/nL))
  We[1:nD,1:nL] <-W 
  for (d in 1:nD){
    for ( p in 1:((nRE/nL)-1)){
      col <- min(which(apply(We,2,sum)==0))
      row <- min(which(apply(We,1,sum)==0))
      w <- W[d,which(mapping_vec==d)]
      We[row,col:(col+length(w)-1)] <- w
    }
   
    
  }
  B_Lambda <- We %*% B_Omega %*% t(We)
  print(B_Lambda)
  A_Omega <- matrix(vec_alpha_ij, nrow = nL, ncol = nL, byrow = TRUE)
  A_Lambda <- W %*% A_Omega %*% t(W)   
  vec_alpha_ij <- as.numeric(t(A_Lambda))  # row-wise
  
  para_weights <- unlist(weights_list)
  
  paraB <- ifelse(pos[5]==pos[6],numeric(0),opt[(pos[5]+1):pos[6]])
  paraSig <- opt[(pos[6]+1):pos[7]]
  ParaTransformY <- opt[(pos[8]+1):pos[9]]
  para_surv <- ifelse(pos[9]==pos[10],numeric(0),opt[(pos[9]+1):pos[10]])
  
  ## ---- (9) Construct return vector replacing Ω-scale with Λ-scale --
  b_inv <- c(
    alpha_mu0,
    alpha_mu,
    as.numeric( B_Lambda[ lower.tri(B_Lambda, diag=TRUE) ] ),
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv,
    para_weights
  )
  
  
  return(b_inv)
}





