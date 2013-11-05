##' Returns liklihood of models relative to the best (based on AIC)
##' 
##' @description Calculates relative liklihood of models based on Akaike Information Criterion
##' @param AIC_df data frame of AIC results, as provided by the AIC function
##' @return The input data frame, with a `liklihood` column

AIC_lik <- function(AIC_df) {
  
  AIC_df$rel_liklihood <- exp((min(AIC_df$AIC)-AIC_df$AIC)/2)
  AIC_df
  
}