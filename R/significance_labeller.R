##' Create 'significance code' based on p values
##' 
##' @Description Returns a nicely-formatted p value (if print_p=TRUE) or a signficance code (e.g. S or NS).
##' @export

significance_labeller <- function(p, alpha=0.05, sig_text="S", nonsig_text="NS", print_p=TRUE, sig_figs=2, tiny=0.001) {
  
  if(print_p) {
    # Create a label of the form "p=x"
    if (p<tiny) { 
      label_text <- paste("p<", tiny, sep="")
    } else {
      label_text <- paste("p=", signif(p, sig_figs), sep="")
    }
    
    
  } else {
    # Create a significance code (e.g. "S", "NS")
    if (p <= alpha) {
      label_text <- sig_text
    } else {
      label_text <- nonsig_text
    }
  }
  
  label_text
}