#' Repeat LP-specific columns in a data frame
#'
#' Creates a new data frame by repeating groups of columns corresponding to LP identifiers
#' (e.g., "DeltaLP.1.", "DeltaLP.2.") a specified number of times for each LP.
#'
#' @param df A \code{data.frame} containing columns named with LP patterns (e.g., "DeltaLP.<id>.<var>").
#' @param repeats_lp A named integer vector where colnames are LP IDs (as character) and values are the number of repetitions.
#'        LPs not listed will default to 1 repetition.
#' @param add_suffix Logical; if \code{TRUE}, adds suffixes like "__rep1", "__rep2" to repeated columns for uniqueness.
#' @param pattern_fun Optional function to extract LP IDs from column colnames. Defaults to regex for "DeltaLP.<id>.".
#'
#' @return A \code{data.frame} with repeated columns according to \code{repeats_lp}.
#'
#' @details
#' Columns are grouped by LP ID based on their colnames. Non-LP columns are kept once by default.
#' The function preserves the original column order and repeats LP blocks contiguously.
#'
#' @examples
#' df <- data.frame(
#'   `DeltaLP.1.(Intercept)` = rnorm(4),
#'   `DeltaLP.1.Time`        = rnorm(4),
#'   `DeltaLP.2.(Intercept)` = rnorm(4),
#'   `DeltaLP.2.Time`        = rnorm(4),
#'   check.colnames = FALSE
#' )
#'
#' # Repeat LP 1 three times, LP 2 twice
#' repeats_lp <- c("1" = 3, "2" = 2)
#' df_new <- repeat_lp_columns(df, repeats_lp, add_suffix = TRUE)
#' colnames(df_new)
#'
#' @export

model_matrix_nL <- function(df, repeats_lp, add_suffix = TRUE,mod=1) {
  

  nms <- colnames(df)
  # Extract LP ids per column (character, NA if no match)
  if(mod==1){
    matches <- as.numeric(sub("^.*?\\.(\\d+).*", "\\1", nms))
  }else{
    matches <-as.numeric(gsub("[^0-9]", "", nms))
  }
 
  
  unique_lp <- unique(matches)
  # Prepare output columns
  out_list <- list()
  
  # Preserve the original column order by iterating over df columns and handling by LP id blocks
 
  for (i in unique_lp) {
    out_list[[i]] <-df[,rep(which(matches==i), times =length(which(repeats_lp==i)))]
    # colnames(out_list[[i]]) <- rep(paste0(nms[which(matches==i)],1:length(which(repeats_lp==i))),each=length(which(repeats_lp==i)))
  }
  
 
  out <- do.call("cbind",out_list)
  return(out)
}
