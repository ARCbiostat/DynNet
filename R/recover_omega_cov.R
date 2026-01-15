#' Recover a single Omega-level covariance matrix from a single alternating Lambda-level matrix
#'
#' @param Blambda  (nD*nTypes) x (nD*nTypes) Lambda-level var/cov (non-Cholesky).
#'                 Rows/cols alternate by type: int0, int, slope (or int0, int).
#' @param mappingL1L2  nD x nL; (d,l) is the weight ω_{d,l} of Omega l in Lambda d.
#' @param method  "structured" (default) or "pseudoinverse".
#' @param rho_int0  scalar or length nD vector: intra-Ω correlation for INT0 within Λ_d group.
#' @param rho_int   scalar or length nD vector: intra-Ω correlation for INT within Λ_d group.
#' @param rho_slope scalar or length nD vector: intra-Ω correlation for SLOPE within Λ_d group.
#'                  Ignored if slopes absent.
#' @param intercept_var_fixed if TRUE, force Var(Λ_d intercept)=1 in the INT block.
#'
#' @return A single Omega-level covariance matrix BOmega of size (nL + 2*nL) x (nL + 2*nL)
#'         ordered as: [all Ω int0] then [Ω1(int), Ω1(slope), Ω2(int), Ω2(slope), ...].
#'         Attributes include a list of components for inspection.
#' @export
recover_omega_cov <- function(Blambda,
                              Blambda_fixed,
                              mappingL1L2,
                              method = c("structured", "pseudoinverse"),
                              rho_int0 = 0,
                              rho_int  = 0,
                              rho_slope = 0) {
  
  method <- match.arg(method)
  
  nB <- nrow(Blambda)
  if (nB != ncol(Blambda)) stop("Blambda must be square.")
  
  # ---- helpers
  pinv <- function(A, tol = .Machine$double.eps) {
    sv <- svd(A); d <- sv$d
    Dplus <- diag(ifelse(d > tol, 1/d, 0), nrow = length(d))
    sv$v %*% Dplus %*% t(sv$u)
  }
  nearest_psd <- function(S, eps = 1e-10) {
    S <- (S + t(S))/2
    ed <- eigen(S, symmetric = TRUE)
    Spsd <- ed$vectors %*% diag(pmax(ed$values, eps)) %*% t(ed$vectors)
    (Spsd + t(Spsd))/2
  }
  
  rn <- rownames(Blambda); cn <- colnames(Blambda)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn))
    warning("Row/col names of Blambda differ; proceeding but consistent names are recommended.")
  mappingL1L2 <- t(mappingL1L2)
  nD <- nrow(mappingL1L2)
  nL <- ncol(mappingL1L2)
  W  <- mappingL1L2
  # --- name-based parsing (preferred)
  idx_int0 <- 1:nD
  idx_int <- setdiff( which(
    grepl("\\(Intercept\\)|Intercept(?!0)", rn, ignore.case = TRUE, perl = TRUE)
  ),idx_int0)
  idx_slope <- setdiff(seq_len(nB), union(idx_int0, idx_int))
  
  has_slopes <- length(idx_slope) > 0

  omega_usage <- colSums(abs(W) > 0)
  overlapped <- which(omega_usage > 1)
  
  # --- extract Λ blocks by type
  BLambda_int0  <- Blambda[idx_int0, idx_int0, drop = FALSE]
  BLambda_int   <- Blambda[idx_int,  idx_int,  drop = FALSE]
  BLambda_slope <- if(length(idx_slope)>0) Blambda[idx_slope, idx_slope, drop = FALSE] else NULL
  BLambda_int_slope <- if (length(idx_slope)>0) Blambda[idx_int, idx_slope, drop = FALSE] else NULL
  
  
  BLambda_fixed_int0  <- Blambda_fixed[idx_int0, idx_int0, drop = FALSE]
  BLambda_fixed_int   <- Blambda_fixed[idx_int,  idx_int,  drop = FALSE]
  BLambda_fixed_slope <- if(length(idx_slope)>0) Blambda_fixed[idx_slope, idx_slope, drop = FALSE] else NULL
  BLambda_fixed_int_slope <- if (length(idx_slope)>0) Blambda_fixed[idx_int, idx_slope, drop = FALSE] else NULL
  
  # broadcast rhos
  if (length(rho_int0) == 1L) rho_int0 <- rep(rho_int0, nD)
  if (length(rho_int)  == 1L) rho_int  <- rep(rho_int,  nD)
  if (!is.null(BLambda_slope) && length(rho_slope) == 1L) rho_slope <- rep(rho_slope, nD)
  
  # --- compute BOmega_* per type
  if (method == "structured") {
    if (length(overlapped) > 0) {
      stop("Structured method assumes disjoint Omega groups per Lambda. Overlaps: ",
           paste(overlapped, collapse = ", "),
           ". Use method='pseudoinverse' or redesign mapping.")
    }
    
    groups <- lapply(seq_len(nD), function(d) which(abs(W[d, ]) > 0))
    
    BOmega_fixed_int0  <- map_lambda_constraint_to_omega(BLambda_fixed_int0,mappingL1L2 = mappingL1L2)
    BOmega_fixed_int  <- map_lambda_constraint_to_omega(BLambda_fixed_int,mappingL1L2 = mappingL1L2)
    BOmega_fixed_slope <- if (!is.null(BLambda_slope)) map_lambda_constraint_to_omega(BLambda_fixed_slope,mappingL1L2 = mappingL1L2) else NULL
    
    # cross-type int–slope (diagonal per Ω, group-specific value)
    BOmega_fixed_int_slope <- if (!is.null(BLambda_slope)) map_lambda_constraint_to_omega(BLambda_fixed_int_slope,mappingL1L2 = mappingL1L2) else NULL
    
    
    BOmega_int0  <- matrix(0, nL, nL)
    BOmega_int   <- matrix(0, nL, nL)
    BOmega_slope <- if (!is.null(BLambda_slope)) matrix(0, nL, nL) else NULL
    # cross-type int–slope (diagonal per Ω, group-specific value)
    BOmega_int_slope <- if (!is.null(BLambda_slope)) matrix(0, nL, nL) else NULL
    
    # INT0
    for (d in seq_len(nD)) {
      G  <- groups[[d]]
      wd <- W[d, G]
      varLambda_d <- BLambda_int0[d, d]
      sum_w2   <- sum(wd^2)
      sum_pair <- ((sum(wd))^2 - sum_w2) / 2
      denom <- sum_w2 + 2 * rho_int0[d] * sum_pair
      if (denom <= 0) stop(sprintf("Non-positive denominator for Lambda %d (int0).", d))
      sigma2_d <- varLambda_d / denom
      block <- matrix(rho_int0[d] * sigma2_d, length(G), length(G))
      diag(block) <- sigma2_d
      BOmega_int0[G, G] <- block
    }
    
    for (d in seq_len(nD)) {
      G  <- groups[[d]]
      wd <- W[d, G]
      varLambda_d <- BLambda_int[d, d]
      sum_w2   <- sum(wd^2)
      sum_pair <- ((sum(wd))^2 - sum_w2) / 2
      denom <- sum_w2 + 2 * rho_int[d] * sum_pair
      if (denom <= 0) stop(sprintf("Non-positive denominator for Lambda %d (int).", d))
      sigma2_d <- varLambda_d / denom
      block <- matrix(rho_int[d] * sigma2_d, length(G), length(G))
      diag(block) <- sigma2_d
      BOmega_int[G, G] <- block
    }
    
    # SLOPE + cross-type INT–SLOPE
    if (!is.null(BLambda_slope)) {
      for (d in seq_len(nD)) {
        G  <- groups[[d]]
        wd <- W[d, G]
        # SLOPE variances/covariances
        varLambda_d <- BLambda_slope[d, d]
        sum_w2   <- sum(wd^2)
        sum_pair <- ((sum(wd))^2 - sum_w2) / 2
        denom <- sum_w2 + 2 * rho_slope[d] * sum_pair
        if (denom <= 0) stop(sprintf("Non-positive denominator for Lambda %d (slope).", d))
        sigma2_d <- varLambda_d / denom
        block <- matrix(rho_slope[d] * sigma2_d, length(G), length(G))
        diag(block) <- sigma2_d
        BOmega_slope[G, G] <- block
        
        # INT–SLOPE covariance (allowed only "within Lambda" and only same Ω -> diagonal per Ω)
        # Match Lambda-level cross covariance: Cov(Λ_int, Λ_slope) = sum_l w_l^2 * tau_d
        covLambda_d <- BLambda_int_slope[d, d]
        tau_d <- if (sum_w2 > 0) covLambda_d / sum_w2 else 0
        BOmega_int_slope[cbind(G, G)] <- tau_d  # place on the diagonal for Ω in group G
      }
    }
    
  } else {
    # pseudoinverse (allows overlaps; we still enforce int–slope "same Ω only" afterwards)
    BLambda_int_eff <- BLambda_int
    if (intercept_var_fixed) diag(BLambda_int_eff) <- 1
    Wplus <- pinv(W)
    BOmega_int0  <- nearest_psd(Wplus %*% BLambda_int0   %*% t(Wplus))
    BOmega_int   <- nearest_psd(Wplus %*% BLambda_int_eff %*% t(Wplus))
    BOmega_slope <- if (!is.null(BLambda_slope))
      nearest_psd(Wplus %*% BLambda_slope %*% t(Wplus)) else NULL
    
    # cross-type int–slope via pseudoinverse then keep only diag by Ω (to satisfy "only within Lambda")
    BOmega_int_slope <- if (!is.null(BLambda_slope)) Wplus %*% BLambda_int_slope %*% t(Wplus) else NULL
    if (!is.null(BOmega_int_slope)) {
      # zero all off-diagonals (keep only Ω-diagonal)
      BOmega_int_slope <- diag(diag(BOmega_int_slope), nrow = nL, ncol = nL)
    }
  }
  
  # --- assemble final BOmega with requested order:
  # 1) all Ω int0 (nL positions)
  # 2) Ω1(int), Ω1(slope), Ω2(int), Ω2(slope), ...
  nTypes_trailing <- if (is.null(BOmega_slope)) 1L else 2L  # trailing types: int (and slope if present)
  pOmega <- nL + nL * nTypes_trailing
  BOmega <- matrix(0, pOmega, pOmega)
  
  # indices in output
  idx_out_int0 <- seq_len(nL)
  idx_out_trail_start <- nL + 1L
  
  idx_out_int   <- integer(nL)
  idx_out_slope <- if (!is.null(BOmega_slope)) integer(nL) else NULL
  for (l in seq_len(nL)) {
    base <- idx_out_trail_start + (l - 1L) * nTypes_trailing
    idx_out_int[l] <- base
    if (!is.null(BOmega_slope)) idx_out_slope[l] <- base + 1L
  }
  
  # place same-type blocks
  BOmega[idx_out_int0, idx_out_int0] <- BOmega_int0
  BOmega[idx_out_int,  idx_out_int]  <- BOmega_int
  if (!is.null(BOmega_slope)) {
    BOmega[idx_out_slope, idx_out_slope] <- BOmega_slope
    # place cross-type int–slope (same Ω only -> diagonal by Ω)
    # This is a diagonal matrix in Ω-space, so we map diagonals to (int,slope) pairs
    for (l in seq_len(nL)) {
      tau_ll <- BOmega_int_slope[l, l]
      if (tau_ll != 0) {
        BOmega[idx_out_int[l], idx_out_slope[l]] <- tau_ll
        BOmega[idx_out_slope[l], idx_out_int[l]] <- tau_ll
      }
    }
  }
  
  #costraints
 
  BOmega_fixed <- matrix(1, pOmega, pOmega)
  
  # indices in output
  idx_out_int0 <- seq_len(nL)
  idx_out_trail_start <- nL + 1L
  
  idx_out_int   <- integer(nL)
  idx_out_slope <- if (!is.null(BOmega_fixed_slope)) integer(nL) else NULL
  for (l in seq_len(nL)) {
    base <- idx_out_trail_start + (l - 1L) * nTypes_trailing
    idx_out_int[l] <- base
    if (!is.null(BOmega_fixed_slope)) idx_out_slope[l] <- base + 1L
  }
  
  # place same-type blocks
  BOmega_fixed[idx_out_int0, idx_out_int0] <- BOmega_fixed_int0
  BOmega_fixed[idx_out_int,  idx_out_int]  <- BOmega_fixed_int
  if (!is.null(BOmega_fixed_slope)) {
    BOmega_fixed[idx_out_slope, idx_out_slope] <- BOmega_fixed_slope
    # place cross-type int–slope (same Ω only -> diagonal by Ω)
    # This is a diagonal matrix in Ω-space, so we map diagonals to (int,slope) pairs
    for (l in seq_len(nL)) {
      tau_ll <- BOmega_fixed_int_slope[l, l]
        BOmega_fixed[idx_out_int[l], idx_out_slope[l]] <- tau_ll
        BOmega_fixed[idx_out_slope[l], idx_out_int[l]] <- tau_ll
      
    }
  }
  
  
  
  omega_names <- paste0("Omega", seq_len(nL))
  
  # Build labels following the output order:
  # [all int0], then [Ω1(int), Ω1(slope), Ω2(int), Ω2(slope), ...]  (or just int if no slopes)
  labels_int0 <- paste0(omega_names, "_int0")
  if (has_slopes) {
    labels_trailing <- as.vector(rbind(
      paste0(omega_names, "_int"),
      paste0(omega_names, "_slope")
    ))
  } else {
    labels_trailing <- paste0(omega_names, "_int")
  }
  labels <- c(labels_int0, labels_trailing)
  
  dimnames(BOmega) <- list(labels, labels)
  dimnames(BOmega_fixed) <- list(labels, labels)
  
  # diagnostics
  offdiag_near_zero <- function(B) {
    if (is.null(B)) return(TRUE)
    A <- B; diag(A) <- 0; all(abs(A) <= 1e-8)
  }
  diagnostics <- list(
    overlapping_omegas = overlapped,
    BLambda_int0_offdiag_near_zero  = offdiag_near_zero(BLambda_int0),
    BLambda_int_offdiag_near_zero   = offdiag_near_zero(BLambda_int),
    BLambda_slope_offdiag_near_zero = offdiag_near_zero(BLambda_slope),
    BLambda_int_slope_offdiag_near_zero = offdiag_near_zero(BLambda_int_slope)
  )
  
  # attach useful components
  attr(BOmega, "components") <- list(
    BOmega_int0      = BOmega_int0,
    BOmega_int       = BOmega_int,
    BOmega_slope     = BOmega_slope,
    BOmega_int_slope = if (exists("BOmega_int_slope")) BOmega_int_slope else NULL,
    W                = W,
    BLambda_int0     = BLambda_int0,
    BLambda_int      = BLambda_int,
    BLambda_slope    = BLambda_slope,
    BLambda_int_slope= BLambda_int_slope,
    diagnostics      = diagnostics,
    method           = method
  )
  
  return(list(BOmega,BOmega_fixed))
 
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
    if (length(G) > 0) {
      for (d2 in seq_len(nD)) {
        if (LambdaConstraint[d, d2] == 1) {
          G2 <- which(abs(mappingL1L2[d2, ]) > 0)
          OmegaConstraint[G, G2] <- 1
        }
      }
    }
  }
  OmegaConstraint
}
