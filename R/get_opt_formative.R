#' Title
#'
#' @param temp 
#'
#' @returns
#' @export
#'
#' @examples
get_opt_formative <- function(temp,paras_block_dim,mapping,mapping2){
  opt <- temp$b
  pos <- c(1,cumsum(paras_block_dim))
  mappingLP2LP1 <- pmin(table(mapping, mapping2), 1)
  nL <- ncol(mappingLP2LP1)
  nD <- nrow(mappingLP2LP1)
  mappingLP2LP1_vec <- apply(mappingLP2LP1,2,function(x)which(x!=0))
  
  alpha_mu0_trans <- opt[pos[1]:pos[2]]
  n_mu0 <- length(alpha_mu0_trans)/nL
  map_mu0 <- rep(c(1:n_mu0),times=nL)
  mappingLP2LP1_vec_mu0 <- rep(mappingLP2LP1_vec,each=n_mu0) 
  alpha_mu0 <- as.numeric(do.call(rbind,lapply(1:n_mu0,function(x)tapply(alpha_mu0_trans[map_mu0==x],mappingLP2LP1_vec_mu0[map_mu0==x],sum))))
  
  alpha_mu_trans <- opt[(pos[2]+1):pos[3]]
  n_mu <- length(alpha_mu_trans)/nL
  map_mu <- rep(c(1:n_mu),times=nL)
  mappingLP2LP1_vec_mu <- rep(mappingLP2LP1_vec,each=n_mu) 
  alpha_mu <- as.numeric(do.call(rbind,lapply(1:n_mu,function(x)tapply(alpha_mu_trans[map_mu==x],mappingLP2LP1_vec_mu[map_mu==x],sum))))
  
  alpha_D <- opt[(pos[3]+1):pos[4]]
  vec_alpha_ij <- opt[(pos[4]+1):pos[5]]
  
  
  
  b_inv <- c(
    alpha_mu0,
    alpha_mu,
    alpha_D,
    vec_alpha_ij,
    paraB,
    paraSig,
    ParaTransformY,
    para_weights,
    para_surv
  )
  
  temp$b <- b_inv
  return(temp)
}