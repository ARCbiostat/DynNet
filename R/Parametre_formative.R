#' Function to initialize parameters k in multivariate DynNet model
#'
#' @param K number of the markers
#' @param nD number of the latent processes
#' @param nL number of exogeneous latent processes
#' @param mapping.to.LP indicates which outcome measured which latent process, it is a mapping table between
#'  outcomes and latents processes
#' @param mapping.to.LP indicates which outcome measured which exogenous latent process, it is a mapping table between
#'  outcomes and enogenous latents processes (only used for structural model)
#' @param vec_ncol_x0n vector of number of columns of model.matrix for baseline's submodel
#' @param n_col_x number of overall columns of model.matrix for change's submodel
#' @param nb_RE number of random effects
#' @param stochErr indicates if the structural model contain stochastique components
#' @param indexparaFixeUser position of parameters to be constrained
#' @param paraFixeUser values associated to the index of parameters to be constrained
#' @param L number of columns of model.matrix for temporal influences model
#' @param paras.ini initial values for parameters, default values is NULL
#' @param ncolMod.MatrixY vector of number of columns of model.matrix for transformation submodel
#' @param link indicates link used to transform outcome
#' @param npara_k marker-specific number of transformation parameters
#' @param Survdata dataset for survival model
#' @param basehaz type of baseline hazard function
#' @param knots_surv knots for splines in baseline hazard (to develop)
#' @param assoc specification of association between longitudinal and survival models
#' @param truncation boolean for delayed entry
#' @param data dataset for longitudinal model
#' @param outcomes names of the outcomes
#' @param df degree of freedom for link==splines
#' @param nE number of survival events
#' @param np_surv cause-specific number of fixed effects in survival models
#' @param fixed.survival.models specification of survival models (without interactions)
#' @param interactionY.survival.models specification of interactions in survival models
#' @param nYsurv number of fixed effects in survival model
#' @return a list
#'
#' @importFrom stats qunif median
#' @importFrom mstate trans.comprisk msprep expand.covs
#' @importFrom survival coxph survreg
#'
#'
Parametre_formative <- function(K,
                       nD,
                       nL,
                       mapping.to.LP,
                       mapping.to.LP2,
                       vec_ncol_x0n,
                       n_col_x,
                       nb_RE,
                       stochErr = FALSE,
                       indexparaFixeUser = NULL,
                       paraFixeUser = NULL,
                       L = 1,
                       paras.ini,
                       ncolMod.MatrixY,
                       link,
                       npara_k,
                       Survdata = NULL,
                       basehaz = NULL,
                       knots_surv = NULL,
                       assoc = NULL,
                       truncation = F,
                       data,
                       outcomes,
                       df,
                       nE = 0,
                       np_surv = 0,
                       fixed.survival.models = NULL,
                       interactionY.survival.models = NULL,
                       nYsurv = 0,
                       names_x,
                       names_x0,
                       names_z,
                       names_z0,
                       names_y,
                       varcovRE.format="cholesky") {
  cl <- match.call()
  
  #   require(MASS)
  #initialisation des parametres
  # L = number of parameters for each coefficient a of matrix A
  # K = number of outcomes
  #======================================================================================
  
  if(varcovRE.format=="cholesky"){
    nb_paraD = nb_RE*(nb_RE+1)/2
  }
  
  if(varcovRE.format=="block"){
    nq <- (nb_RE-nD)/nD
    # random int + cholesky + rho int +rho int slopes
    nb_paraD=nD+(nq*(nq+1)/2)*nD+(nD^2-nD)/2+nD*nq
    
  }
  indexparaFixeForIden <- NULL
  # if user not specified initial parameters
  
  if (is.null(basehaz))
    basehaz <- "Weibull" #not to have NULL value in C++ code
  
  if (is.null(paras.ini)) {
    stop("Initial parameters need to specified!")
  }
  
  if(!is.list(paras.ini)){
    stop("Formative structural model requires initial parameters object to be initialized using enter_param()")
  }
  
  # if user specified initial parameters

    map_p <-list() #attribution of parameters to latent process
    
    p <- 0 # position in the initialize parameters
    cpt1 <-0 # counter for parameterd
    cpt2<-0 # loop counter
    
    #alpha_mu0
    alpha_mu0 <- paras.ini$alpha_mu0
    map_p$alpha_mu0 <- as.integer(sub("^LP0\\.(\\d+).*", "\\1", names_x0))
    names(map_p$alpha_mu0) <- paste0("alpha_mu0", 1:length(alpha_mu0))
    
    #counting
    p <- p+ sum(vec_ncol_x0n)
    index_paraFixe_mu0_constraint <-NULL
    for(n in 1:nD){
      #alpha_mu0[(cpt2+1)] <- 0
      cpt2 <- cpt2 + vec_ncol_x0n[n]
      cpt1 <- cpt1 + vec_ncol_x0n[n]
    }
    paraFixe_mu0_constraint <- rep(1,nD)
    
    #alpha_mu
    alpha_mu <- paras.ini$alpha_mu
    map_p$alpha_mu <- as.integer(sub("^DeltaLP\\.(\\d+).*", "\\1", names_x))
    names(map_p$alpha_mu) <- paste0("alpha_mu", 1:length(alpha_mu))
    
    p <- p+n_col_x
    cpt1 <- cpt1 + n_col_x
    
    #random effects
    alpha_D <- paras.ini$alpha_D
    if(varcovRE.format=="cholesky"){
      alpha_D_matrix <- DparChol(nb_RE,alpha_D)
     
    }
    
    if(varcovRE.format=="block"){
      alpha_D_matrix <-DparBlock(nD,(nb_RE-nD)/nD,alpha_D)
      print(alpha_D_matrix)
    }
    
    
    RE_z0 <- as.integer(sub("\\(.*\\)", "", names_z0))
    RE_z <- as.integer(sub("\\(.*\\)", "", names_z))
    colnames(alpha_D_matrix) <- c(names_z0,names_z)
    rownames(alpha_D_matrix) <- c(names_z0,names_z)
    
    #counting
    to_nrow <- nb_RE
    i_alpha_D <- 0
    index_paraFixeDconstraint <- NULL
    
    for(n in 1:nD){
      #if(link[n] != "thresholds")
      #alpha_D[i_alpha_D+1] <- 1
      i_alpha_D <- i_alpha_D + to_nrow
      cpt1 <- cpt1 + to_nrow
      to_nrow <- to_nrow -1
    }
    p <- p+nb_paraD
    paraFixeDconstraint <- rep(1,nD)
   
    # para of transition matrix vec_alpha_ij
    vec_alpha_ij <- paras.ini$vec_alpha_ij
    map_p$vec_alpha_ij <- rep(1:nD, each = nD)# to check if correct
    names(map_p$vec_alpha_ij) <- paste0("vec_alpha_ij", 1:length(vec_alpha_ij))
    
    p <- p + L*nD*nD
    cpt1 <- cpt1 + L*nD*nD
    
    # paraB
    paraB <- NULL
    if (stochErr == TRUE) {
      paraB <- paras.ini$paraB
      map_p$paraB <- 1:nD
      names(map_p$paraB) <- paste0("paraB", 1:length(paraB))
      
      p <- p + nD
      cpt1 <- cpt1 + nD
    }
    
    #paraSig
    paraSig <- paras.ini$paraSig
    map_p$paraSig <- mapping.to.LP
    names(map_p$paraSig)<- paste0("paraSig", 1:length(paraSig))
    
    p <- p + K
    cpt1 <- cpt1 + K
    
    ### para of link function
    ParaTransformY <- paras.ini$ParaTransformY
    i_para <- 0
    for (k in 1:K) {
      if (link[k] == "linear" & ParaTransformY[i_para + 2] == 0) {
        stop('Second parameter for linear link function cannot be set at 0 (variance)')
      }
      i_para <- i_para + npara_k[k]
    }
    
    
    map_p$ParaTransformY <- rep(mapping.to.LP, times = as.numeric(table(sub(
      "\\..*", "", names_y
    ))))
    names(map_p$ParaTransformY) <- paste0("ParaTransformY", 1:length(ParaTransformY))
    
    cpt1 <- cpt1 + ncolMod.MatrixY
    p <- p + ncolMod.MatrixY
    
    
  #Survival
  knots_surv <- c(0,0)
  para_surv <- paras.ini$para_surv
  
  if(!is.null(Survdata)){
   
    np_baz <- ifelse(basehaz=="Weibull",2, 0)# changer 0!!
    
    for (jj in 1:nE){
     
      p <- p + np_baz  # change here?
  
      p <- p + np_surv[jj] # change here?
    }
    if(basehaz=="Splines") cat('add number of parameters for splines in p and para_surv')
    if(basehaz=="Splines") cat('Define knots_surv para_basehaz')
    
  }
  
  
  W_raw <- matrix(0, nrow=nD, ncol=nL)
  for(d in 1:nD){
    idx <- which(mappingLP2LP1[,d] == 1)
    k <- length(idx)
    if(k > 0){
      w <- runif(k, 0.1, 1)   # enforce imbalance
      w <- w / sum(w)
      W_raw[d, idx] <- w
    }
  }
  
  print("initial weights")
  print(W_raw)
  
  para_weights <- apply(W_raw,2,sum)
  #final vector of initial parameters
  paras <- c(
    alpha_mu0,
    alpha_mu,
    alpha_D,
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv,
    para_weights
  )
  
  if(!is.null(paras.ini)){
    if( length(paras) != p ){
      stop("The length of paras.ini is not correct.")
    }}
  
  
  paras_length <- unlist(lapply(list(
    alpha_mu0,
    alpha_mu,
    alpha_D,
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_surv
  ), length))
  return(
    list(
      para = paras,
      paraOpt = paraOpt,
      paraFixe = paraFixe,
      posfix = posfix,
      L = L,
      basehaz = basehaz,
      knots_surv = knots_surv,
      np_surv = np_surv,
      assoc = assoc,
      truncation = truncation,
      nb_paraD= nb_paraD,
      paras_length=paras_length
    )
  )
}
