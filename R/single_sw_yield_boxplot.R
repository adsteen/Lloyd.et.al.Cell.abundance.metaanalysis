##' Make a single yield boxplot (seds and sw, despite the name)
##' 
##' @description Function to make generic boxplots as in the supplemental figures
##' @export

single_sw_yield_boxplot <- function(d, xvar, yvar, sig_text=NA, jitter_width=0.25, alpha=0.25, ps=0.75, text_size=4) {
  
  text_df <- data.frame(x=-Inf, y=Inf, label=sig_text)
  
  p <- ggplot(d, aes_string(x=xvar, y=yvar)) + 
    geom_point(position=position_jitter(jitter_width), alpha=alpha, size=ps) +
    geom_boxplot(alpha=0, outlier.size=NA) +
    geom_text(data=text_df, aes(x=x, y=y, label=label), vjust=1.2, hjust=-0.2, size=text_size) +
    aes(ymin=0) +
    theme(axis.text.x=element_text(angle=-45, hjust=0),
          plot.margin=unit(c(5, 5, 5, 5), "mm"),
          panel.margin=unit(c(5, 5, 5, 5), "mm")) 
  p
}