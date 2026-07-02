#' Recover a single Omega-level covariance matrix from a single alternating Lambda-level matrix
#'
#' @param Blambda  (nD*nTypes) x (nD*nTypes) Lambda-level var/cov (non-Cholesky).
#'                 Rows/cols alternate by type: int0, int, slope (or int0, int).
#' @param mappingL1L2  nD x nL; (d,l) is the weight ω_{d,l} of Omega l in Lambda d.
#' @return A single Omega-level covariance matrix BOmega of size (nL + 2*nL) x (nL + 2*nL)
#'         ordered as: [all Ω int0] then [Ω1(int), Ω1(slope), Ω2(int), Ω2(slope), ...].
#'         Attributes include a list of components for inspection.
#' @export
recover_omega_cov <- function(alpha_D,
                              nb_RE,
                              mappingL1L2) {
  
  
  # nB <- nrow(Blambda)
  # if (nB != ncol(Blambda)) stop("Blambda must be square.")
  # 
  # rn <- rownames(Blambda); cn <- colnames(Blambda)
  # if (!is.null(rn) && !is.null(cn) && !identical(rn, cn))
  #   warning("Row/col names of Blambda differ; proceeding but consistent names are recommended.")
  # 
  nD <- nrow(mappingL1L2)
  nL <- ncol(mappingL1L2)
  W  <- matrix(mappingL1L2,ncol=nL,nrow=nD)
  
  
  alpha_D[1:nD] <- 
    apply(W, 1, function(w){
      sum(w^2) +
        2 * 0 * sum(outer(w, w)[upper.tri(outer(w, w))])
    })
  Blamdabyp <- DparBlock(nD,(nb_RE-nD)/nD,alpha_D)
  has_slope <- ifelse(nrow(Blamdabyp)>nD*2,TRUE,FALSE)
  Blambda   = reorder_bytype(Blamdabyp ,nD=nD,has_slope=has_slope)
  
  # # --- name-based parsing (preferred)
  # idx_int0 <- 1:nD
  # idx_int <- setdiff( which(
  #   grepl("\\(Intercept\\)|Intercept(?!0)", rn, ignore.case = TRUE, perl = TRUE)
  # ),idx_int0)
  # idx_slope <- setdiff(seq_len(nB), union(idx_int0, idx_int))
  # 
  # has_slope <- length(idx_slope) > 0

  omega_usage <- colSums(abs(W) > 0)
  overlapped <- which(omega_usage > 1)

    if (length(overlapped) > 0) {
      stop("Omega groups per Lambda needs to be disjoint")
    }
    
    if(has_slope)A <- as.matrix(Matrix::bdiag(W, W, W)) else A <- as.matrix(Matrix::bdiag(W, W))
    Ap <- MASS::ginv(as.matrix(A))
    
    BOmega <- Ap %*% Blambda %*% t(Ap)
    BOmega <- (BOmega + t(BOmega))/2
    eps <- 1e-3
    BOmega <- BOmega + eps * diag(nrow(BOmega))
    sd0 <- sqrt(diag(BOmega))
    R0  <- cov2cor(BOmega)
    alpha <- 0.7   # o 0.5
    
    R1 <- R0
    R1[row(R1) != col(R1)] <- alpha * R1[row(R1) != col(R1)]
    BOmega <- diag(sd0) %*% R1 %*% diag(sd0)
    
    diag(BOmega)[1:nL] <- 1
    print("eerore ricostruzione")
    print(max(abs(
      as.numeric(A %*% BOmega %*% t(A)) -
        as.numeric(Blambda)
    ))>0.01)
    if (
      max(abs(
        as.numeric(A %*% BOmega %*% t(A)) -
        as.numeric(Blambda)
      ))>0.01
    ){
      Delta <- Blambda - A %*% BOmega %*% t(A)
      print(
        max(abs(Delta))
        
      )
      print(range(abs(Delta)))
      print(Delta)
      # stop("something wrong mapping BOmega")
    }
      
    print("check invertibility")
    
    print(eigen(BOmega)$values)
  # attach useful components
  attr(BOmega, "components") <- list(
    W                = A,
    Wp=Ap
  )
  
  return(list(BOmega,1:nL))
 
}



# Block-wise Cholesky 
chol_by_block <- function(BOmega, blocks, adjust_2x2_if_needed = TRUE) {
  stopifnot(is.matrix(BOmega), nrow(BOmega) == ncol(BOmega))
  p <- nrow(BOmega)
  if (length(blocks) != p) stop("Length of 'blocks' must equal dimension of BOmega.")
  
  # Group indices by block ID
  block_groups <- split(seq_len(p), blocks)
  
  # Initialize result
  R <- matrix(0, nrow = p, ncol = p)
  
  # Process each block
  for (bname in names(block_groups)) {
    idx <- block_groups[[bname]]
    Bsub <- BOmega[idx, idx, drop = FALSE]
    k <- length(idx)
    
    if (k == 1L) {
      # Scalar case
      v <- Bsub[1, 1]
      if (v <= 0) v <- .Machine$double.eps
      R[idx, idx] <- sqrt(v)
    } else if (k == 2L) {
      # 2x2 case with optional PD correction
      det2 <- Bsub[1,1]*Bsub[2,2] - Bsub[1,2]^2
      if (adjust_2x2_if_needed && (Bsub[1,1] <= 0 || Bsub[2,2] <= 0 || det2 <= 0)) {
        v1 <- max(Bsub[1,1], .Machine$double.eps)
        v2 <- max(Bsub[2,2], .Machine$double.eps)
        tau_max <- sqrt(v1*v2) * (1 - 1e-8)
        tau <- max(min(Bsub[1,2], tau_max), -tau_max)
        Bsub <- matrix(c(v1, tau, tau, v2), 2, 2)
      }
      Rsub <- tryCatch(chol(Bsub), error = function(e) NULL)
      if (is.null(Rsub)) stop(sprintf("Block '%s' not PD even after adjustment.", bname))
      R[idx, idx] <- Rsub
    } else {
      # Larger block
      Rsub <- tryCatch(chol(Bsub), error = function(e) NULL)
      if (is.null(Rsub)) stop(sprintf("Block '%s' not PD.", bname))
      R[idx, idx] <- Rsub
    }
  }
  
  # to get the lower triangular matrix
  t(R)
}



# Create Omega-level constraint matrix from Lambda-level constraint matrix
# Inputs:
#   LambdaConstraint: nD x nD matrix with 1 where parameter is fixed
#   mappingL1L2: nD x nL mapping
# Output:
#   OmegaConstraint: nL x nL matrix with 1 where corresponding Omega pairs are fixed
map_lambda_constraint_to_omega <- function(LambdaConstraint, mappingL1L2) {
  
  nD <- nrow(mappingL1L2)
  nL <- ncol(mappingL1L2)
  OmegaConstraint <- matrix(0, nL, nL)
  
  for (d in seq_len(nD)) {
    G <- which(abs(mappingL1L2[d, ]) > 0)
    
    for (d2 in seq_len(nD)) {
      
      if (LambdaConstraint[d, d2] == 1) {
        G2 <- which(abs(mappingL1L2[d2, ]) > 0)
        
        if (d == d2) {
         
          OmegaConstraint[cbind(G, G)] <- 1
          
        } else {
          
          common <- intersect(G, G2)
          if (length(common) > 0) {
            OmegaConstraint[cbind(common, common)] <- 1
          }
        }
      }
    }
  }
  
  OmegaConstraint
}


reorder_bytype <- function(B, nD, has_slope = TRUE) {
  stopifnot(is.matrix(B), nrow(B) == ncol(B))
  
  nTypes <- if (has_slope) 3 else 2
  if (nrow(B) != nD * nTypes) {
    stop("Dimensions do not match nD * nTypes")
  }
  
  # --- rebuild forward permutation EXACTLY
  idx_int0 <- 1:nD
  idx_int  <- (nD + 1):(2 * nD)
  
  if (has_slope) {
    idx_slope <- (2 * nD + 1):(3 * nD)
    perm_trailing <- as.vector(rbind(idx_int, idx_slope))
  } else {
    perm_trailing <- idx_int
  }
  
  perm_forward <- c(idx_int0, perm_trailing)
  
 
  perm_inv <- order(perm_forward)
  
  B_new <- B[perm_inv, perm_inv]
  B_new
}

reorder_byprocess <- function(B, nD, has_slope = TRUE) {
  stopifnot(is.matrix(B), nrow(B) == ncol(B))
  
  nTypes <- if (has_slope) 3 else 2
  if (nrow(B) != nD * nTypes) {
    stop("Dimensions do not match nD * nTypes")
  }
  
  idx_int0 <- 1:nD
  idx_int  <- (nD + 1):(2 * nD)
  
  if (has_slope) {
    idx_slope <- (2 * nD + 1):(3 * nD)
    perm_trailing <- as.vector(rbind(idx_int, idx_slope))
  } else {
    perm_trailing <- idx_int
  }
  
  perm <- c(idx_int0, perm_trailing)
  
  B_new <- B[perm, perm]
  attr(B_new, "perm") <- perm  # optional
  B_new
}