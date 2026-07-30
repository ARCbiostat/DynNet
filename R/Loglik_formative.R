Loglik_formative <- function(K , nD, mapping,nL,mapping2,paraOpt,  paraFixe , posfix , paras_k ,
                             sequence , type_int , ind_seq_i, MCnr , nmes ,
                             m_is , Mod_MatrixY , Mod_MatrixYprim , df,
                             x ,x_omega, z , q , nb_paraD ,
                             x0 ,x0_omega, z0 , q0 , varcov_format ,
                             data_surv , data_surv_intY , nYsurv , basehaz , knots_surv, 
                             np_surv , survival , assoc , truncation, 
                             nE, Xsurv1 , Xsurv2 ,
                             if_link , zitr, ide,
                             tau , tau_is, 
                             modA_mat, DeltaT, ii,paras_dim,paras_length,
                             map_p,
                             varcovRE.format){
  
  mappingLP2LP1 <- pmin(table(mapping2, mapping), 1)
  
  para <- rep(NA,sum(length(paraOpt),length(paraFixe)))
  para[posfix==0] <- paraOpt
  para[posfix==1] <- paraFixe
  
  para_ends <- cumsum(paras_length)
  para_starts <- para_ends - paras_length + 1
  # For types with zero elements
  para_starts[paras_length == 0] <- NA
  para_ends[paras_length == 0] <- NA
 
  # here I start transforming parameters
  
  # Weigths for formative structural model
  para_weights_raw <- para[para_starts["para_weights"]:para_ends["para_weights"]]
  
  W_raw <- matrix(0, nrow=nD, ncol=nL)
  for(d in 1:nD){
    idx <- which(mappingLP2LP1[,d] == 1)
    k <- length(idx)
    if(k > 0){
      eta <- para_weights_raw[idx]
      w <- exp(eta)
      w <- w / sum(w)
      
      W_raw[d,idx] <- w
    }
  }
  
  
  nb_RE <- sum(c(q,q0))
  
  mappingLP2LP1_weights_raw <- t(mappingLP2LP1) * W_raw
  
  alpha_mu0 <- para[para_starts["alpha_mu0"]:para_ends["alpha_mu0"]]
  alpha_mu <- para[para_starts["alpha_mu"]:para_ends["alpha_mu"]]
  alpha_D <- para[para_starts["alpha_D"]:para_ends["alpha_D"]]
  vec_alpha_ij <- para[para_starts["vec_alpha_ij"]:para_ends["vec_alpha_ij"]]
  if(!is.na(para_starts["paraB"])) paraB <- para[para_starts["paraB"]:para_ends["paraB"]]else paraB <- NULL
  paraSig <- para[para_starts["paraSig"]:para_ends["paraSig"]]
  ParaTransformY <- para[para_starts["ParaTransformY"]:para_ends["ParaTransformY"]]
  if(!is.na(para_starts["para_surv"]))para_surv <- para[para_starts["para_surv"]:para_ends["para_surv"]]else para_surv <- NULL
  
  
  indexFixe_alpha_mu0 <- which(posfix[para_starts["alpha_mu0"]:para_ends["alpha_mu0"]]==1)
  indexFixe_alpha_mu <- which(posfix[para_starts["alpha_mu"]:para_ends["alpha_mu"]]==1)
  indexFixe_alpha_D <- which(posfix[para_starts["alpha_D"]:para_ends["alpha_D"]]==1)
  indexFixe_vec_alpha_ij <- which(posfix[para_starts["vec_alpha_ij"]:para_ends["vec_alpha_ij"]]==1)
  if(!is.na(para_starts["paraB"])) indexFixe_paraB <- qwhich(posfix[para_starts["paraB"]:para_ends["paraB"]]==1)else indexFixe_paraB <- NULL
  indexFixe_paraSig <- which(posfix[para_starts["paraSig"]:para_ends["paraSig"]]==1)
  indexFixe_ParaTransformY <- which(posfix[para_starts["ParaTransformY"]:para_ends["ParaTransformY"]]==1)
  if(!is.na(para_starts["para_surv"]))indexFixe_para_surv <- which(posfix[para_starts["para_surv"]:para_ends["para_surv"]]==1)else indexFixe_para_surv <- NULL
  
  alpha_D_matrix <- DparBlock(nL,(nb_RE-nL)/nL,alpha_D)
  
  
  if(varcovRE.format=="cholesky"){
    stop("Formative model and cholesky not developed.")
    indexFixe_alpha_D_matrix <- DparChol(nb_RE,indexFixe_alpha_D)
    
    
    #alpha_D
    alpha_D_matrix_trans <- recover_omega_cov(
      Blambda   = alpha_D_matrix,
      Blambda_fixed = indexFixe_alpha_D_matrix,
      mappingL1L2   = mappingLP2LP1_weights
    )
    
    reps <- length(row.names(alpha_D_matrix_trans[[1]])[-c(1:nL)])/nL
    blocks <- c(1:nL,rep(apply(mappingLP2LP1,1,function(x)which(x!=0)),times=reps)+(nL+1))
    
    
    nb_paraD <- nrow(alpha_D_matrix_trans[[1]])*(nrow(alpha_D_matrix_trans[[1]])+1)/2
    
    
    
    #chol_D_block <- chol_by_block(alpha_D_matrix_trans[[1]],blocks = blocks) I can't do it because it is not consistent with C++
    # print("original")
    # print(alpha_D_matrix_trans[[1]])
    chol_D <- t(chol(alpha_D_matrix_trans[[1]]))
    # print("chol")
    # print(chol_D)
    alpha_D_trans <- as.numeric(chol_D[lower.tri(chol_D, diag = TRUE)])
    # print("numeric version")
    # print(alpha_D_trans)
    #check
    # print("check transformation")
    # print(identical(DparChol(nrow(alpha_D_matrix_trans[[1]]),alpha_D_trans),alpha_D_matrix_trans[[1]]))
    # print(DparChol(nrow(alpha_D_matrix_trans[[1]]),alpha_D_trans))
    # print(alpha_D_matrix_trans[[1]])
    
    
    
    indexFixe_alpha_D_trans_matrix <- alpha_D_matrix_trans[[2]]
    indexFixe_alpha_D_matrix <- indexFixe_alpha_D_trans_matrix[lower.tri(indexFixe_alpha_D_trans_matrix, diag = TRUE)]
    indexFixe_alpha_D_trans <- which(indexFixe_alpha_D_matrix==1)
    
  }
  
  
  
  # if(varcovRE.format=="block"){
  #   
  #   #indexFixe_alpha_D_matrix <- DparBlock(nD,(nb_RE-nD)/nD,indexFixe_alpha_D) # Not used
  #   has_slope <- ifelse(nrow(alpha_D_matrix)>nD*2,TRUE,FALSE)
  #   if(nrow(alpha_D_matrix)>(nD*3)){
  #     stop("Random effects with a non-linear effect of time are currently not implemented")
  #   }
  # 
  #   
  #   
  #   alpha_D_matrix_trans <- recover_omega_cov(
  #     alpha_D = alpha_D,
  #     nb_RE = nb_RE,
  #     mappingL1L2   = mappingLP2LP1_weights_raw
  #   )
  #   
  #   
  #   print(alpha_D_matrix_trans[[1]])
  #   nq <- (nrow(alpha_D_matrix_trans[[1]])-nL)/nL
  #   # random int + cholesky + rho int +rho int slopes
  #   nb_paraD=nL+(nq*(nq+1)/2)*nL+(nL^2-nL)/2+nL*nq
  #   print("number par random effects")
  #   print(nb_paraD)
  #   
  #   # change order by process
  #   alpha_D_matrix_trans_byp <- reorder_byprocess(alpha_D_matrix_trans[[1]],nL,has_slope=has_slope)
  #   
  #   var.int <- diag(alpha_D_matrix_trans[[1]])[1:nL]
  #   varcovRE.time <-lapply(seq_len(nL), function(d) {
  #     
  #     if (has_slope) {
  #       base <- nL + 2*(d - 1)
  #       idx <- c(base + 1, base + 2)   # int_d, slope_d
  #     } else {
  #       idx <- nL + d                  # int_d only
  #     }
  #     
  #     alpha_D_matrix_trans_byp[idx, idx, drop = FALSE]
  #   })
  #   
  #   rho.int <- sapply(1:round(nL/2), function(l) {
  #     idx_int0 <- l
  #     idx_int <- (l+1):nL
  #     
  #     cov <- alpha_D_matrix_trans_byp[idx_int0, idx_int]
  #     var0 <- diag(alpha_D_matrix_trans_byp)[idx_int0]
  #     var1 <- diag(alpha_D_matrix_trans_byp)[idx_int]
  #     
  #     
  #     
  #     rho <- cov / sqrt(var0 * var1)
  #     
  #     
  #   },simplify = T)
  #   
  #   print("rho.int")
  #   print(rho.int)
  #   
  #   rho.int.time <-  sapply(seq_len(nL), function(l) {
  #     
  #     idx_int0 <- l
  #     
  #     idx_int <- if (has_slope) {
  #       nL + 2*(l - 1) + 1
  #     } else {
  #       nL + l
  #     }
  #     
  #     cov <- alpha_D_matrix_trans_byp[idx_int0, idx_int]
  #     var0 <- alpha_D_matrix_trans_byp[idx_int0, idx_int0]
  #     var1 <- alpha_D_matrix_trans_byp[idx_int,  idx_int]
  #     
  #     cov / sqrt(var0 * var1)
  #   })
  #   
  #   print("rho.int.time")
  #   print(rho.int.time)
  #   varcovRE.time.chol <- lapply(varcovRE.time, function(x) t(chol(x))[lower.tri(t(chol(x)),diag = T)])
  #   alpha_D_trans <- c(var.int,unlist(varcovRE.time.chol),inv_rho(unlist(rho.int)),unlist(lapply(rho.int.time,function(x)inv_rho(x))))
  #   print(alpha_D_trans)
  #   print(DparBlock(nL,nq,alpha_D_trans))
  #   
  #   # indexFixe_alpha_D_trans_matrix_byp <- reorder_byprocess(alpha_D_matrix_trans[[2]],nL,has_slope=has_slope)
  #   # 
  #   # var.int.fix <- diag(alpha_D_matrix_trans[[2]][1:nL,1:nL])
  #   # varcovRE.time.fix <-lapply(seq_len(nL), function(d) {
  #   # 
  #   #   if (has_slope) {
  #   #     base <- nL + 2*(d - 1)
  #   #     idx <- c(base + 1, base + 2)   # int_d, slope_d
  #   #   } else {
  #   #     idx <- nL + d                  # int_d only
  #   #   }
  #   # 
  #   #   indexFixe_alpha_D_trans_matrix_byp[idx, idx, drop = FALSE]
  #   # })
  #   # 
  #   # rho.int.fix <- sapply(seq_len(nL), function(l) {
  #   # 
  #   #   idx_int0 <- l
  #   # 
  #   #   idx_int <- if (has_slope) {
  #   #     nL + 2*(l - 1) + 1
  #   #   } else {
  #   #     nL + l
  #   #   }
  #   # 
  #   #   indexFixe_alpha_D_trans_matrix_byp[idx_int0, idx_int]
  #   # 
  #   # })
  #   # 
  #   # rho.int.time.fix <-  sapply(seq_len(nL), function(l) {
  #   # 
  #   #   idx_int0 <- l
  #   # 
  #   #   idx_int <- if (has_slope) {
  #   #     nL + 2*(l - 1) + 1
  #   #   } else {
  #   #     nL + l
  #   #   }
  #   # 
  #   #   indexFixe_alpha_D_trans_matrix_byp[idx_int0, idx_int]
  #   # 
  #   # })
  #   
  #   
  #   
  #   # indexFixe_alpha_D_trans <-c(var.int.fix,
  #   #                             unlist(lapply(varcovRE.time.fix, function(x) as.numeric(x[lower.tri(x,diag = T)]))),
  #   #                             rho.int.fix,
  #   #                             rho.int.time.fix)
  #   
  #   indexFixe_alpha_D_trans <- 1:nL
  # }
  
  #alrady on scale of Omegas
  indexFixe_alpha_D_trans <- indexFixe_alpha_D
  alpha_D_trans <- alpha_D
  alpha_D_matrix_trans <- alpha_D_matrix
  # ===============================
  ## Scaling 
  ## ===============================
  R0 <- alpha_D_matrix_trans[1:nL, 1:nL]
  
  B_lambda00 <- W_raw %*% R0 %*% t(W_raw)
  
  D <- diag(
    1 / sqrt(diag(B_lambda00))
  )
  
  W_scaled <- D %*% W_raw
  mappingLP2LP1_weights <- t(mappingLP2LP1) * W_scaled
  para_weights <-  apply(W_scaled,2,sum)
  
  #alpha_mu0
  alpha_mu0_trans <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_p$weights==x],function(l)alpha_mu0[map_p$alpha_mu0==x]*l)))
  names(alpha_mu0_trans) <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_p$weights==x],function(l)names(map_p$alpha_mu0)[map_p$alpha_mu0==x])))
  
  if(length(indexFixe_alpha_mu0)!=0){
    mapping_par <- match(names(alpha_mu0_trans), names(map_p$alpha_mu0))
    
    indexFixe_alpha_mu0_trans <- unlist(lapply(indexFixe_alpha_mu0,function(x)which(mapping_par==x)))
    paraFixe_alpha_mu0_trans <- alpha_mu0_trans[indexFixe_alpha_mu0_trans]
  }else{
    indexFixe_alpha_mu0_trans <- NULL
    paraFixe_alpha_mu0_trans <- NULL
  }
  
  #alpha_mu
  
  alpha_mu_trans <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_p$weights==x],function(l)alpha_mu[map_p$alpha_mu==x]*l)))
  names(alpha_mu_trans) <- unlist(lapply(1:nD,function(x)lapply(para_weights[map_p$weights==x],function(l)names(map_p$alpha_mu)[map_p$alpha_mu==x])))
  
  if(length(indexFixe_alpha_mu)!=0){
    mapping_par <- match(names(alpha_mu_trans), names(map_p$alpha_mu))
    
    indexFixe_alpha_mu_trans <- unlist(lapply( alpha_mu,function(x)which(mapping_par==x)))
    paraFixe_alpha_mu_trans <- alpha_mu0_trans[indexFixe_alpha_mu_trans]
  }else{
    indexFixe_alpha_mu_trans <- NULL
    paraFixe_alpha_mu_trans <- NULL
  }
  
  
  #vec_alpha_ij
  vec_alpha_ij_trans <-as.vector(t(invert_vec_alpha_ij(matrix(vec_alpha_ij,nrow=nD,byrow = T),t(mappingLP2LP1_weights))))
  indexFixe_vec_alpha_ij_trans <- if(length(indexFixe_vec_alpha_ij)>0)map_fixed_lambda_to_omega(indexFixe_vec_alpha_ij,mappingLP2LP1)else integer(0)
  
  paras_trans <- c(
    alpha_mu0_trans,
    alpha_mu_trans,
    alpha_D_trans,
    vec_alpha_ij_trans,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv
  )
  
  print(paras_trans)
  
  paras_trans_length <- unlist(lapply(list(
    alpha_mu0_trans,
    alpha_mu_trans,
    alpha_D_trans,
    vec_alpha_ij_trans,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv,
    para_weights
  ), length))
  
  indexparaFixeUser_trans <- c(indexFixe_alpha_mu0_trans,
                               indexFixe_alpha_mu_trans,
                               indexFixe_alpha_D_trans,
                               indexFixe_vec_alpha_ij_trans,
                               indexFixe_paraB,
                               indexFixe_paraSig,
                               indexFixe_ParaTransformY,
                               indexFixe_para_surv)
  
  
  paraFixeUser_trans <- c(alpha_mu0_trans[indexFixe_alpha_mu0_trans],
                          alpha_mu_trans[indexFixe_alpha_mu_trans],
                          alpha_D_trans[indexFixe_alpha_D_trans],
                          vec_alpha_ij_trans[indexFixe_vec_alpha_ij_trans],
                          paraB[indexFixe_paraB],
                          paraSig[indexFixe_paraB],
                          ParaTransformY[indexFixe_ParaTransformY],
                          para_surv[indexFixe_para_surv])
  
  
  indexparaFixeUser_trans <- c(indexFixe_alpha_mu0_trans,
                               indexFixe_alpha_mu_trans+paras_trans_length[1],
                               indexFixe_alpha_D_trans+sum(paras_trans_length[1:2]),
                               indexFixe_vec_alpha_ij_trans+sum(paras_trans_length[1:3]),
                               indexFixe_paraB+sum(paras_trans_length[1:4]),
                               indexFixe_paraSig+sum(paras_trans_length[1:5]),
                               indexFixe_ParaTransformY+sum(paras_trans_length[1:6]),
                               indexFixe_para_surv+sum(paras_trans_length[1:7]))
  
  
  #initialisation
  #   paraOpt <- paras
  posfix <- rep(0, length(paras_trans)) # 0 = non fixe 1 = fixe # initialisation
  # constraining of parameters==============
  indexFixe <- NULL # not used
  
  if (!is.null(indexparaFixeUser_trans)) {
    
    indexFixe <- sort(unique(c(indexFixe, indexparaFixeUser_trans)))
  }
  paraFixe <- rep(NA, length(posfix))
  if (!is.null(paraFixeUser_trans)) {
    paraFixe[c(indexparaFixeUser_trans)] <- paraFixeUser_trans
  }
  
  #not used
  # paraFixe[index_paraFixe_mu0_constraint] <- rep(0, K)
  # paraFixe[index_paraFixeDconstraint] <- rep(1, K)
  if (sum(!is.na(paraFixe)) == 0) {
    paraFixe = -1
    paraOpt <- paras_trans
  } else{
    paraFixe <- paraFixe[!is.na(paraFixe)]
    posfix[indexFixe] <- 1 # fixation des paras d'indexes dans indexparaFixe
    paras_trans[indexFixe] <- paraFixe
    paraOpt <- paras_trans[-indexFixe]
  }
  
  # print("Final check parametre formative")
  # print("paras")
  # print(paras_trans)
  # print(length(paras_trans))
  # 
  print("paraOpt")
  print(paraOpt)
  print(length(paraOpt))
  
  print("paraFixe")
  print(paraFixe)
  print(length(paraFixe))
  
  
  
 
    if(type_int == 1){
      sequence  <- randtoolbox::sobol(n = MCnr, dim = sum(length(q0)+length(q)), scrambling = 1, normal = TRUE, init=T)
    }else if(type_int == 2){
      sequence  <- randtoolbox::halton(n = MCnr, dim = sum(length(q0)+length(q)), normal = TRUE, init=T) 
    }else if(type_int == 3){
      sequence  <- randtoolbox::torus(n = MCnr, dim = sum(length(q0)+length(q)), normal = TRUE, init=T) 
    }else{
    sequence  <- matrix(0,nL,1)
    }
  
  check_arma_compatibility <- function(...) {
    args <- list(...)
    for (name in names(args)) {
      obj <- args[[name]]
      cat(name, ": ")
      
      if (is.null(obj)) {
        cat("NULL\n")
      } else if ((is.numeric(obj) || is.logical(obj)) && length(obj) == 1) {
        if (is.logical(obj)) {
          cat("scalar bool candidate (OK)\n")
        } else {
          cat("scalar numeric candidate (OK)\n")
        }
      } else if (is.numeric(obj) && is.vector(obj)) {
        cat("arma::vec candidate, length =", length(obj), "\n")
      } else if (is.numeric(obj) && is.matrix(obj)) {
        cat("arma::mat candidate, dim =", paste(dim(obj), collapse = " x "), "\n")
      } else {
        cat("WARNING: unexpected type =", class(obj), "\n")
      }
    }
  }
  
  # # fix model.matrix for formative
  # mappingLP2LP1 <- pmin(table(mapping, mapping2), 1)
  # mappingLP2LP1_vec <- apply(mappingLP2LP1,2,function(x)which(x!=0))
  # x <-model_matrix_nL(x,repeats_lp = mappingLP2LP1_vec)
  # x0 <-model_matrix_nL(x0,repeats_lp = mappingLP2LP1_vec)
  
  
  # LogLik
  res <- Loglik(
    K = K,
    nD = nL,
    mapping =  as.numeric(mapping2),
    paraOpt = as.numeric(paraOpt),
    paraFixe = paraFixe,
    posfix = posfix,
    paras_k = paras_k,
    sequence = sequence,
    type_int = type_int,
    ind_seq_i = ind_seq_i,
    MCnr = MCnr,
    nmes = nmes,
    m_is = m_is,
    Mod_MatrixY = Mod_MatrixY,
    Mod_MatrixYprim = Mod_MatrixYprim,
    df = df,
    x = as.matrix(x_omega),
    z = as.matrix(z),
    q = q,
    nb_paraD = nb_paraD,
    x0 = as.matrix(x0_omega),
    z0 = as.matrix(z0),
    q0 = q0,
    varcov_format=varcov_format,
    data_surv = data_surv,
    data_surv_intY = data_surv_intY,
    nYsurv = nYsurv,
    basehaz = basehaz,
    knots_surv = knots_surv,
    np_surv = np_surv,
    survival = survival,
    assoc =  assoc,
    truncation = truncation,
    nE = nE,
    Xsurv1 = Xsurv1,
    Xsurv2 = Xsurv2,
    if_link = if_link,
    zitr = zitr,
    ide = ide,
    tau = tau,
    tau_is = tau_is,
    modA_mat = modA_mat,
    DeltaT=DeltaT,
    ii = ii
  )
  
  return(res)
  
}