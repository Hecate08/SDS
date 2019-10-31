#' Score c
#'
#' This function calculates the clinical independence score based on attribute values and clinical data.
#' @param attrib Matrix of atrribute values with attribute names as row names and sample names as column names.
#' @param clin data.frame of clinical data with dummy variables for categorical variables and with sample names as row names and variable names as column names.
#' @keywords feature selection
#' @return Vector of scores with attribute names as names. Higher is better.
#' @export
#' @examples
#' get_score_c()
get_score_c <- function(attrib,clin){
  lenE <- dim(attrib)[1]
  R2adj <- rep(NA,lenE)
  for(f in 1:lenE){
    fit <- lm(attrib[f,] ~ .,data = clin)
    R2adj[f] <- summary(fit)$adj.r.squared
  }
  R2adj_ok <- as.numeric(gsub(NaN, 1, R2adj))
  names(R2adj_ok) <- rownames(attrib)

  score2 <- 1 - R2adj_ok

  return(score2)
}
