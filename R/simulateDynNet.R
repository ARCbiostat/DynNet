
simulateDynNet <- function(Ni,
                           structural.model,
                           measurement.model, 
                           parameters, 
                           TempsFin,
                           DeltaT,
                           basehaz = NULL, 
                           seed=NULL, 
                           varcovRE.format="cholesky",...){
  
  
  
  if(!missing(seed))
    set.seed(seed)
  
  ### check if all component of the model specification are well filled ####
  if(missing(structural.model))stop("The argument structural.model must be specified")
  if(missing(measurement.model))stop("The argument measurement.model must be specified")
  if(missing(parameters))stop("The argument parameters must be specified")
  
  
  if(is.null(structural.model$fixed.DeltaLP))stop("The argument structural.model$fixed.DeltaLP must be specified")
  #if(is.null(structural.model$random.DeltaLP))stop("The argument structural.model$random.DeltaLP must be specified")
  if(is.null(structural.model$trans.matrix))stop("The argument structural.model$trans.matrix must be specified")
  if(is.null(structural.model$delta.time)){
    structural.model$delta.time <- 1
  }
  
  if(is.null(measurement.model$link.functions) || all(is.null(measurement.model$link.functions$links))){
    links <- NULL
    knots <- NULL
    measurement.model$link.functions =list(links = links, knots = knots)
  }else{
    if(all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="quant", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="manual", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="equi", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="linear", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="thresholds", x))) ==0))
      stop("The only available link functions are 'linear', 'splines' and 'thresholds' functions.")
  }
  if(is.null(parameters$Fixed.para.index))stop("The argument parameters$Fixed.para.index cannot be NULL")
  if(is.null(parameters$Fixed.para.values))stop("The argument parameters$Fixed.para.values cannot be NULL")
  
  if(is.null(option$makepred)){
    option$makepred <- F
  }
  
  
  survival= FALSE
  if(!is.null(structural.model$fixed.survival)){
    survival = TRUE
    fixed.survival <- structural.model$fixed.survival
  }
  
  ### identification of model components #####
  ## components of structural model
  fixed_X0 <- structural.model$fixed.LP0
  fixed_DeltaX <- structural.model$fixed.DeltaLP
  randoms_DeltaX <- structural.model$random.DeltaLP
  mod_trans <- structural.model$trans.matrix
  DeltaT <- structural.model$delta.time
  
  
  interactionY.survival <- NULL
  if(!is.null(structural.model$interactionY.survival)){
    interactionY.survival <- structural.model$interactionY.survival
  }
  
  # components of measurement model
  link <- measurement.model$link.functions$links
  knots <- measurement.model$link.functions$knots
  ## components of parameters initialisation
  
  if(is.null(attr(parameters,"components")))warning("Object of initial values was not created with the enter_param() function. Note that this is recommanded.")
  if(!is.null(attr(parameters,"components")) & !identical(attr(parameters,"components")$structural.model,structural.model)) stop("The object of initial values was created with a different structural model specification.")
  if(!is.null(attr(parameters,"components")) & !identical(attr(parameters,"components")$measurement.model,measurement.model)) stop("The object of initial values was created with a different measurement model specification.")
  indexparaFixeUser <- parameters$Fixed.para.index
  paraFixeUser <- parameters$Fixed.para.values
  paras.ini <- parameters$paras.ini
  
  if(!inherits(fixed_DeltaX,"formula")) stop("The argument fixed_DeltaX must be a formula")
  
  ### outcomes and latent processes ####
  outcome <- as.character(attr(terms(fixed_DeltaX),"variables"))[2]
  outcomes_by_LP<-strsplit(outcome,"[|]")[[1]]
  nD <- length(outcomes_by_LP) # nD: number of latent process
  
  formative <- any(grepl("\\([^)]*[+][^)]*\\)", outcomes_by_LP))
  
 
  outcomes <- NULL
  mapping.to.LP <- NULL
  mapping.to.LP2 <- NULL
  for(n in 1:nD){
    outcomes_n <- strsplit(outcomes_by_LP[n],"[+]")[[1]]
    outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
    outcomes_n <- unique(outcomes_n)
    if(is.null(outcomes_n)) stop("at least one marker must be specified for each latent process" )
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  
  if(formative){
    outcomes <- as.character(sapply(outcomes,FUN = function(x)gsub("[()+]","",x),simplify = FALSE))
    for (n in 1:nD){
      outcomes_n <-  regmatches(outcomes_by_LP[n], gregexpr("\\([^()]*\\)|Y\\d+", outcomes_by_LP[n]))[[1]]
      outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[()]","",x),simplify = FALSE))
      outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
      outcomes_n <- unique(outcomes_n)
      for (l in 1:length(outcomes_n)){
        outcomes_n2 <- regmatches(outcomes_n[l], gregexpr("Y\\d+", outcomes_n[l]))[[1]]
        if(is.null(outcomes_n2)) stop("at least one marker must be specified for each latent exogeneous process" )
        add <- ifelse(n==1,0,max(mapping.to.LP2))
        mapping.to.LP2 <- c(mapping.to.LP2, rep(l+add,length(outcomes_n2)))
      }
      
    }
  }
  outcomes <- sub("[()]", "", outcomes)
  
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  if (formative) {
    nL <- max(mapping.to.LP2)
  } else {
    nL <- NULL
  }
  
  
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  if(nD !=length(fixed_DeltaX.models)) stop("The number of models does not correspond to the number of latent processes")
  
  if(formative){
    if(nL > K){
      stop("There are too many latent processes compared to the indicated number of markers")
    }
  }else{
    if(nD > K){
      stop("There are too many latent processes compared to the indicated number of markers")
    }
  }
  
  ### pre-traitement of fixed effect on initial levels of processes
  if(is.null(fixed_X0)){
    fixed_X0<- ~1
    fixed_X0.models <- rep("1",nD)
  }
  if(!inherits(fixed_X0,"formula")) stop("The argument fixed_X0 must be a formula")
  
  fixed_X0.models =strsplit(gsub("[[:space:]]","",as.character(fixed_X0)),"~")[[2]]
  fixed_X0.models<- as.vector(strsplit(fixed_X0.models,"[|]")[[1]]) 
  for(nd in 1:nD){
    if(fixed_X0.models[nd]=="~-1") fixed_X0.models[nd] <-"~1" # au moins l'intcpt
  }
  
  ### pre-traitement of fixed effect on survival 
  fixed.survival.models <- NULL
  truncation = FALSE
  assoc <- 0
  if(survival){
    if(is.null(fixed.survival)){
      fixed.survival<- ~1
    }else{
      if(!inherits(fixed.survival,"formula")) stop("The argument fixed.survival must be a formula") 
    }
    
    fixed.survival.models <- strsplit(gsub("[[:space:]]","",as.character(fixed.survival)),"~")[[2]]
    covsurv <- unique(as.vector(strsplit(fixed.survival.models,"[|*+]")[[1]]))
    fixed.survival.models <- as.vector(strsplit(fixed.survival.models,"[|]")[[1]]) 
   
    
    if(!is.null(option$assocT)){
      if(!option$assocT%in%c("r.intercept", "r.slope", "r.intercept/slope", "c.value"))
        stop("assocT should be defined as r.intercept, r.slope, r.intercept/slope, c.value.")
      assoc <- switch(option$assocT, "r.intercept"=0, "r.slope"=1, "r.intercept/slope"=2, "c.value"=3, "c.slope"=4, "c.value/slope"=5)
    }else{
      assocT <- "r.intercept/slope"
      assoc <- 2 # random intercept and slope
    }
    
    if(assoc <= 2){
      message(" add interactions ui * X in survival model")
      
    }
    if(!is.null(option$truncation)){
      truncation <- option$truncation
    }
  }
  
  ### pre-traitement of interactions with Y on survival 
  
  interactionY.survival.models <- NULL
  if(survival){
    
    if(!is.null(interactionY.survival)){
      if(!inherits(interactionY.survival,"formula")) stop("The argument interactionY.survival must be a formula") 
      interactionY.survival.models <- strsplit(gsub("[[:space:]]","",as.character(interactionY.survival)),"~")[[2]]
      intYsurv <- (as.vector(strsplit(interactionY.survival.models,"[|*+]")[[1]]))
      interactionY.survival.models <- as.vector(strsplit(interactionY.survival.models,"[|]")[[1]]) 
     
      if(any(grepl( "*", interactionY.survival, fixed = TRUE)))stop("Only + terms should be included in interactionY.survival, no *.")
    }
  }
  
  
  ### pre-traitement of random effect on processes  intercept and slope
  #### randoms effet on DeltaLP 
  randoms_X0.models <- rep("1",nD)
  #### randoms effet on DeltaX
  
  if(missing(randoms_DeltaX) || is.null(randoms_DeltaX)){
    randoms_DeltaX<- ~1
    randoms_DeltaX.models <- rep("1",nD)
  }
  if(!inherits(randoms_DeltaX,"formula")) stop("The argument random must be a formula")
  randoms_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(randoms_DeltaX)),"~")[[2]]
  randoms_DeltaX.models<-strsplit(randoms_DeltaX.model,"[|]")[[1]]    
  
  #### traitement of  mod_trans: transition matrix##### 
  if(missing(mod_trans)){
    mod_trans <- ~ 1 # constant transition matrix
  } 
  if(!inherits(mod_trans,"formula")) stop("The argument mod_trans must be a formula")
  mod_trans.model=strsplit(gsub("[[:space:]]","",as.character(mod_trans)),"~")[[2]]
  
  if(nD!=length(fixed_X0.models)){
    stop("The number of models for initial latent processes does not correspond with the number of latent processes")
  }
  if(nD!=length(fixed_DeltaX.models)){
    stop("The number of models for the change over time of latent processes does not correspond with the number of latent processes")
  }
  
  ### traitement of transformation models ##
  if(is.null(link)){
    link <- rep("linear",K)
  }
  else if(length(link)!=K) stop("The number transformation links must be equal to the number of markers")
  
  
  #################### here starts simulation part #####################
  
  
  # All times and visits assuming cont time observation
  Time <- seq(from = 0, to = TempsFin, by = DeltaT)
  Visit <- seq(0,(length(Time)-1))
  
  ############ data part ############################
  
  subject <- "id"
  Time <- "time"
  data <- data.frame(id=1:Ni,
                     time=Time,
                     visit=Visit)
  I <- length(Visit)
  K <- length(outcomes)
  
      #all predictor of the model==============================================================
     
      all.pred.fixed_X0 <- NULL
      all.pred.fixed_DeltaX <- NULL
      all.pred.randoms_X0 <- NULL
      all.pred.randoms_DeltaX <- NULL
      all.pred.mod_trans <- NULL
      
      for( n in 1: nD){
        all.pred.fixed_X0 <- c(all.pred.fixed_X0,list(strsplit(fixed_X0.models[n],"[+]")[[1]]))
        all.pred.fixed_DeltaX <- c(all.pred.fixed_DeltaX,list(strsplit(fixed_DeltaX.models[n],"[+]")[[1]]))
        all.pred.randoms_X0 <- c(all.pred.randoms_X0,list(strsplit(randoms_X0.models[n],"[+]")[[1]]))
        all.pred.randoms_DeltaX <- c(all.pred.randoms_DeltaX,list(strsplit(randoms_DeltaX.models[n],"[+]")[[1]]))
      }
      #
      all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
      all.pred.fixed_X0<-sapply(all.pred.fixed_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
      #
      all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
      all.pred.fixed_DeltaX<-sapply(all.pred.fixed_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
      #
      all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
      all.pred.randoms_X0<-sapply(all.pred.randoms_X0,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
      #
      all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
      all.pred.randoms_DeltaX<-sapply(all.pred.randoms_DeltaX,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
      #
      all.pred.mod_trans <- c(all.pred.mod_trans,list(strsplit(mod_trans.model,"[+]")[[1]]))
      all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE)
      all.pred.mod_trans<-sapply(all.pred.mod_trans,FUN = function(x)gsub("[*]",":",x),simplify = FALSE)
      #ajout
      all.preds<-unlist(unique(c(unlist(all.pred.fixed_X0), unlist(all.pred.fixed_DeltaX), 
                                 unlist(all.pred.randoms_X0), unlist(all.pred.randoms_DeltaX),
                                 all.pred.mod_trans)))
      all.preds<-inclu_intercerpt(all.preds) # remplace 1 par "(Intercept)"
      all.preds<-unique(unlist(sapply(all.preds,FUN=function(x)strsplit(x,":")[[1]])))
      
      ###all.pred_san_inter signifie all.pred_sans_intercept
      all.preds<-all.preds[-which(all.preds %in% c("(Intercept)"))]
      all.preds <- c(all.preds,Time)
      
      
      cov <- data.table::data.table()
      for (j in unique(c(all.pred.fixed_X0,all.pred.fixed_DeltaX))){
        x[[j]] <- rnorm(Ni)
      }
      
      q0 <- NULL
      nb_paraDw <- 0
      for(n in 1:nD){
        r<-as.formula(paste(subject,randoms_X0.models[n], sep="~-1+"))
        z0n<-model.matrix(r,data=data)
        if(length(z0n)==0){
          col <- paste(n,"zero",sep="")
          z0n<-matrix(assign(col,rep(0,Ni)))
        }
        colnames<-colnames(z0n)
        colnames<-paste(n,colnames,sep="")
        colnames(z0n) <-colnames
        col_n <-c(col_n,list(colnames))
        z0<-cbind(z0,z0n)
        q0 <- c(q0,ncol(z0n))
      }
      
      z <- NULL
      col_n<-list()
      q <- NULL
      nb_paraDu <- 0
      
      for(n in 1:nD){
        r<-as.formula(paste(subject,randoms_DeltaX.models[n], sep="~-1+"))
        zn<-model.matrix(r,data=data)
        if(length(zn)==0){
          col <- paste(n,"zero",sep="")
          zn<-matrix(assign(col,rep(0,Ni)))
        }
        colnames<-colnames(zn)
        colnames<-paste(n,colnames,sep="")
        colnames(zn) <-colnames
        col_n <-c(col_n,list(colnames))
        z<-cbind(z,zn)
        q <- c(q,ncol(zn))
      }
    
  }
  
   #################
  if(varcovRE.format=="block"){
    nb_RE <- length(q0)+length(q)
    matD <- DparBlock(nD,(nb_RE-nD)/nD,unlist(parameters$paraB))
  }
  
  res <- NULL
  nRE <- parameters
  for(i in 1:I){ #cycle over individuals
    
    Y_ij <- NULL
    X_ij_all <-NULL
    
    data_covariate_i <- as.data.frame(cbind(data[i,],cov[i,]))
    cols_data_covariate_i <- colnames(data_covariate_i)
    # model matrix?
    out <- create_x0_x_z0_z_modA_mat(
      data = data_covariate_i,
      fixed_X0.models = fixed_X0.models,
      fixed_DeltaX.models = fixed_DeltaX.models,
      randoms_X0.models = randoms_X0.models,
      randoms_DeltaX.models = randoms_DeltaX.models,
      mod_trans.model = mod_trans.model,
      subject = subject,
      Time = Time
    )
    
    
    #xi pour l'individu i
    x0i <- as.matrix(out$x0)
    xi <- as.matrix(out$x)
    #zi pour l'individu i
    z0i <- as.matrix(out$z0)
    zi <- as.matrix(out$z)
    # modA_mat_i pour l'individu i
    modA_mat_i <- as.matrix(out$modA_mat)
    # we generate the random effects
    Re <- mvrnorm(n = 1, mu = rep(0,(sum(q)+K)), Sigma = matD)
    
    wi <- Re[1:K]
    ui <- Re[(K+1):(sum(q)+K)]
    
    # calcul de X0
    X_ij <- x0i%*%para_mu0 + z0i%*%wi
    X_ij_all <- t(X_ij)
    #measurament errors: allowed K!=ny
    ny = length(M)
    eps <- matrix(NA, length(Time), ny)
    eps[1, ] <- mvrnorm(n = 1, mu = rep(0,ny), Sigma = Sig)
    
    ## boucle pour la recurrence
    X_ij_1 <- X_ij
    
    Aj_1 <- ConstrA(K=K, t=0, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    
  
  
    # latent process
    for(j in 1:(length(Time)-1)){ 
      # il faut sauter xi[1:K,] qui correspond a xi0
      X_ij <- DeltaT*xi[(j*K+1):((j+1)*K),]%*%para_mu + DeltaT*zi[(j*K+1):((j+1)*K),]%*%ui + Aj_1%*% X_ij_1
      eps[j+1,] <- mvrnorm(n = 1, mu = rep(0,ny), Sigma = Sig)
      X_ij_1 <- X_ij
      X_ij_all <- rbind(X_ij_all, t(X_ij))
      Aj_1 <- ConstrA(K=K, t=j, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
    }
    
    # here it changes with formative model
    
    if(formative){
      #ex latent processes
      O_ij_all <- matrix(NA,nrow=nrow( X_ij_all),ncol=max(L))
      for (l in 1:max(L)){
        O_ij_all[,l] <- X_ij_all[,L[l]]*w[l]
      }
    }
    
    
    
    #transformation 
    
    
    y <- NULL
    
    for(m in 1:ny){
      if(links[m]=="linear"){
        y <- cbind(y,(O_ij_all[,L[m]]+eps[,m])*paraEtha1[m] + paraEtha0[m])
        
      } else if (links[m]=="thresholds"){
        
        ## passer des parametres aux thresholds
        thresholds2 <- lapply(thresholds, function(x){c(x[1], x[1]+cumsum(x[-1]^2))})
        
        ##fonction pour passer de H(Y) a Y
        transfinv <- function(ytilde,k)
        {
          nam <- names(thresholds2)
          nam <- tolower(nam)
          namk <- substr(nam[k],1,6)
          
          linktype <- pmatch(namk,c("linear","spline","thresh"))
          
          y <- NA
          
          if(linktype==3) ## thresholds
          {
            v <- c(ytilde,thresholds2[[m]])
            indic <- c(1,rep(0,length(thresholds2[[m]])))
            indic <- indic[order(v)]
            pos <- which(indic==1)
            y <- modalites[[m]][pos]
          }
          
          if(linktype==1) ##linear
          {
            y <- thresholds[[m]][2]*ytilde+thresholds[[m]][1]
          }
          
          if(linktype==2) ##splines
          {
            z <- modalites[[m]] #spline nodes
            ff <- function(x,hy,z,b){transfo_spl(x,z,b)-hy}
            y <- uniroot(f=ff,lower=z[1],upper=z[length(z)],
                         hy=ytilde,z=z,b=thresholds[[m]])$root
          }
          
          
          return(y)
        }
        
        # take the inverse to obtain simulations of Y
        yk <- sapply(O_ij_all[,L[m]]+eps[,m],transfinv,k=m) 
        y <- cbind(y,yk)
      }
    }
    
    Y_ij <- as.data.frame(y)
    colnames(Y_ij)<- paste("Y",1:ny, sep="")
    ### data complet
    data_i_cplt <- as.data.frame(cbind(data_covariate_i,Y_ij))
    data_i_cplt <- data_i_cplt[which(Time%%DeltaTestim==0),]
    data <- rbind(data,data_i_cplt)
    
    
    #survival part
    if(nEvent > 0){
      
      Tevt <- rep(NA,nEvent)
      Devt <- rep(NA,nEvent)
      
      
      xs <- C1[1]
      nweib <- 2
      age0  <- 0
      t.max <- max(Time)#20
      X <- C1[1]
      
      for (j in 1:nEvent){
        
        expSurv <- exp(xs*betaSurv[j])
        
        if(sharedtype == "re"){    # shared random effects
          stop('to develop')
          
          if(!is.null(weibull))
          {
            weibull_j <- weibull[(j-1)*nweib+1:nweib]
            if(nweib==2)
            {
              ## attention prm doivent etre tels que H(t) = w1 * t^w2 (ie parametrisation de logscale=TRUE dans lcmm)
              tempsSurv <- rweibull(1, shape=weibull_j[2], scale=(weibull_j[1]*expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom))^(-1/weibull_j[2]))
              #cat("tempsSurv=",tempsSurv,"\n")
            }
            else
            {
              ## ici prm tels que H(t) = (w1*t)^w2
              unif <- runif(1)
              tempsSurv <- ((-log(unif)/(expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom)))^(1/weibull_j[2]))/weibull_j[1] + weibull_j[3]
            }
          }
          
          if(!is.null(piecewise))
          {
            zl <- piecewise[[paste("nodes",j,sep="")]]
            bl <- piecewise[[paste("brisq",j,sep="")]]
            
            surv <- function(t,zl,bl,expS,p)
            {
              j <- which.max(zl[which(zl<=t)])
              if(j==1) som <- 0 
              else som <- sum(bl[1:(j-1)]*(zl[2:j]-zl[1:(j-1)]))
              
              if(j<length(zl)) surv <- exp(-(som+bl[j]*(t-zl[j]))*expS)
              else surv <- exp(-som*expS)
              
              return(surv-p)
            }
            
            unif <- runif(1)
            zero <- try(uniroot(surv,interval=c(zl[1],zl[length(zl)]),
                                zl=zl,bl=bl,
                                expS=expSurv*exp(t(association[(j-1)*nRE+1:nRE])%*%brandom),
                                p=unif),
                        silent=TRUE)
            if(class(zero)=="try-error") tempsSurv <- Inf
            else tempsSurv <- zero$root
            
          }
        }
        else if(sharedtype == "value"){    # shared latent process current level    #TS
          # code inspired from simulateJM() fct from package JM simulating data with survival model adjusted on the unique outcome current level
          # principle : compute cumulative risk fct values according to several times of event and stop when equals a randomly generated uniform value 
          
          # if(!is.null(piecewise))
          #   stop("piecewise with sharedtype=value not programmed yet")
          
          # if(!is.null(weibull)){
          #   if(nweib != 2) stop("3 parameters Weibull risk not programmed yet")
          weibull_j <- weibull[(j-1)*nweib+1:nweib]
          
          u <- runif(1)  #TS: sim 1 valeur U[0,1] par sujet
          
          
          ###################################################################################
          # function to compute the inverse survival / cumulative risk function
          invS <- function (t, u) {
            
            TD <- function (v) {
              
              ##################
              #Re <- mvrnorm(n = 1, mu = rep(0,(sum(q)+K)), Sigma = matD)
              
              v_ord      <- v[order(v)]
              v_ord_grid <- sapply(v_ord, function(x) round(x/DeltaT))
              Ytild <- matrix(0, length(v_ord_grid),K)
              
              wi <- Re[1:K]
              ui <- Re[(K+1):(sum(q)+K)]
              
              # calcul de X0
              X_ij <- x0i%*%para_mu0 + z0i%*%wi
              X_ij_all <- t(X_ij)
              jj=1
              
              while(v_ord_grid[jj]==0 & jj <= length(v_ord_grid)){
                Ytild[jj,] <- X_ij_all
                jj <- jj+1
              }
              #eps <- matrix(NA, length(Time), K)
              #eps[1, ] <- mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig)
              #Y_ij <- rbind(Y_ij, as.numeric(paraEtha0 + matEtha1%*%(X_ij + eps[1, ] )))
              
              ## boucle pour la recurrence
              X_ij_1 <- X_ij
              
              Aj_1 <- ConstrA(K=K, t=0, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
              
              for(s in 1:(max(v_ord_grid))){ 
                # il faut sauter xi[1:K,] qui correspond a xi0
                X_ij <- DeltaT*xi[(s*K+1):((s+1)*K),]%*%para_mu + DeltaT*zi[(s*K+1):((s+1)*K),]%*%ui + Aj_1%*% X_ij_1
                #eps[j+1,] <- mvrnorm(n = 1, mu = rep(0,K), Sigma = Sig)
                
                #Y_ij <- rbind(Y_ij, as.numeric(paraEtha0 + matEtha1%*%(X_ij + eps[j+1,])))
                #maj X_ij_1 et A_j_1
                X_ij_1 <- X_ij
                X_ij_all <- rbind(X_ij_all, t(X_ij))
                Aj_1 <- ConstrA(K=K, t=s, DeltaT = DeltaT, vec_alpha_ij= vec_alpha_ij, modA_mat=modA_mat_i)
                
                while(v_ord_grid[jj]==s & jj<= length(v_ord_grid)){
                  Ytild[jj,]<-t(X_ij)
                  jj=jj+1
                }
              }
              
              out <- Ytild %*% association[(j-1)*nRE+1:nRE]
              
              return(out)
            }
            
            h <- function (s) {
              # instant risk fct for patient i at time s
              
              TD.i <- TD(s)
              
              #print("Make sure weibull[1] = lambda^alpha and weibull[2] = alpha")
              instant <- weibull_j[1] * weibull_j[2] * s^(weibull_j[2]-1) * exp(TD.i) * expSurv
              
              return(instant)
            }
            
            return(integrate(h, lower = 0, upper = t)$value + log(u))  # cumulative risk fct value + log(uniform value)
          }
          ###################################################################################
          
          risqcum_min =-1 #always negative if age0=0
          #risqcum_min <- invS(age0, u=u) 
          
          risqcum_max <- invS(age0+t.max, u=u)
          
          if(risqcum_min > 0)
          {
            tempsSurv <- age0-1 # event before entry
          }
          else
          {
            if(risqcum_max < 0)
            {
              tempsSurv <- age0+t.max+1 # event after the end of the study
            }
            else
            {
              
              ## find event time between 0 and t.max
              Root <- try(expr = uniroot(invS,
                                         interval = c(1e-05, t.max),
                                         u = u)$root,
                          silent = TRUE)
              
              if(inherits(Root, "try-error"))
              {
                tempsSurv <- as.numeric(2*t.max) #stop("unable to generate a survival time")
              }
              else
              {
                tempsSurv <- Root
              }
            }
          }
          
          
        }
        
        if(tempsSurv < age0) break #on sort de la boucle j
        
        Tevt[j] <- tempsSurv                
        
      }#fin boucle j (1:nEvent)
      
      if(tempsSurv < age0) next # on passe au sujet suivant (boucle i)
      
      tsurvdizaine = F
      if(tsurvdizaine==FALSE)
      {
        Tevt_CR <- min(Tevt,max(Time))
        whichD <- which.min(c(Tevt, max(Time)))
      }
      else
      {
        Tevt_CR <- min(Tevt, max(Time)/10)
        whichD <- which.min(c(Tevt,max(Time)/10))
      }
      
      if(whichD > nEvent){
        D <- 0
      }else{
        D <- whichD
      }
      
      d <- data[which(data$Time<Tevt_CR),] 
      d <- cbind(d,xs, rep(Tevt_CR,dim(d)[1]), rep(D,dim(d)[1]), rep(age0,dim(d)[1]))
      colnames(d) <- c(names(data), "xs","Tevt","Status", "Tentry")
      res <- rbind(res,d)
      
    }else{
      
      res <- rbind(res, data)
    }
  }
  
  

  

