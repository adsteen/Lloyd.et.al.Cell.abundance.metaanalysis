##' AOV permutation test
##' 
##' @description Performs a permutation test to assess variability among groups in data with non-normally-distributed variance
##' @param n Number of permutations
##' @param d Data frame to permute
##' @param cat_var Name of categorical (independant) variable, as a string. Should work as a column number as well.
##' @param data_var Name of random (dependant) variable, as a string. Should work as a column number as well.
##' @param print_hist If `true`, plot the histogram of F values.
##' @param save_hist If `true`, save the histogram of F values
##' @param title Title of the histogram
##' @param fn Filename, as a string, to save the histogram
##' @export

aov_perm_test <- function(n, d, cat_var, data_var, print_hist=FALSE, 
                          save_hist=TRUE, title=NULL, fn=NA) {
  
  if(!is.factor(d[, cat_var])) {
    stop("cat_var needs to refer to a factor")
  }
  if (!is.numeric(d[, data_var])) {
    stop("data_var needs to refer to a numeric column")
  }
  
  y <- d[ , data_var]
  cat <- d[ , cat_var]
  aov_non_perm <- aov(y~cat)
  f_non_perm <- summary(aov_non_perm)[[1]]$"F value"[1]
  
  # initialize the vector of f_values based on permuted data sets
  f_vec <- rep(NA, n)
  
  # Perform the permutation
  for (i in 1:n) {
    aov_perm <- aov(y ~ sample(cat))
    f_vec[i] <- summary(aov_perm)[[1]]$"F value"[1]
  }
  
  # Detgermine limits for a potential plot
  if (f_non_perm > max(f_vec)) {
    xlimits <- c(0, f_non_perm)
  } else {
    xlimits <- c(0, max(f_vec))
  }
  
  if (print_hist) {
    # Plot the histogram
    hist(f_vec, xlim = xlimits, 
         main=title,
         xlab="F value")
    abline(v=f_non_perm, col="red")
    
  }
  
  if (save_hist) {
    if(is.na(fn)) {
      fn <- paste("plots/", Sys.time())
    }
    png(filename=fn, height=3, width=4, units="in", res=300, type="cairo") 
    # Plot the histogram
    hist(f_vec, xlim = xlimits, 
         main=title,
         xlab="F value")
    abline(v=f_non_perm, col="red")
    dev.off()
  }
  
  p.val <- sum(f_vec >= f_non_perm) / length(f_vec)
  #print(paste("f_non_perm is", f_non_perm))
  
  p.val
}