##' Makes plot and calculates summary stats for qPCR using or not using 516 as a primer
##' 
##' @param corrected_seds data frame of corrected seds
##' @return Returns a list with the following components:
##' 
##' * `p_516`, the plot object
##' * med_and_IQR, the data frame of median yield and interquartile range for each depth and primer set (uses or doesn't use 516)
##' * t_test_516 An object of class `htest`, containing the results of the student's t test of yield for samples using vs not using 516 (all depths)
##' * n_points The number of data points for each combination of depth and uses/doesn't use 516

qPCR_516_evaluation <- function(corrected_seds) {
  eval_516 <- corrected_seds[!is.na(corrected_seds$percentqPCR) &
                               !is.na(corrected_seds$Uses.516.for.Arc) &
                               !is.na(corrected_seds$depth_log10), ]
  
  # Calculate median and interquartile range for each method
  med_and_IQR <- ddply(eval_516, .(Uses.516.for.Arc), summarise, 
                       median.percent.arc=median(percentqPCR), 
                       twenty.fifth.percentile=quantile(percentqPCR, probs=0.25),
                       seventy.fifth.percentile=quantile(percentqPCR, probs=0.75))
  
  # Note that a t test isn' really appropriate for this clearly-non-normally distributed data, but whatever. The result is visually obvious anyway
  t_test_516 <- t.test(x=eval_516$percentqPCR[eval_516$Uses.516.for.Arc], y=eval_516$percentqPCR[!eval_516$Uses.516.for.Arc])
  
  # Set qualitative depth bins
  depth_labs <- c("1 cm or less", "1-10 cm", "10 cm-1 m", "1-10 m", "10-100 m", ">100 m")
  eval_516$qual_depth <- cut(eval_516$depth_log10, breaks=-3:3, labels=depth_labs, ordered_result=TRUE)
  
  # Change order of factor levels to make it look better
  eval_516$qual_depth <- factor(eval_516$qual_depth, levels=levels(eval_516$qual_depth)[nlevels(eval_516$qual_depth):1], ordered=TRUE)
  eval_516$Uses.516.for.Arc <- factor(eval_516$Uses.516.for.Arc, levels=c("TRUE", "FALSE"), ordered=TRUE)
  
  # Calculate the number of points for e
  n_points <- ddply(eval_516, c("qual_depth", "Uses.516.for.Arc"), summarise, n.points=length(percentqPCR))
  med_and_IQR <- ddply(eval_516, c("qual_depth", "Uses.516.for.Arc"), summarise, 
                                   median.percentqPCR = median(percentqPCR, na.rm=TRUE),
                                   IQR.percentqPCR = IQR(percentqPCR, na.rm=TRUE))
  
  pdf("~/Dropbox/Metadata Analysis/Revision for AEM/516_evaluation.pdf")
  grid.table(med_and_IQR)
  dev.off()
  
  # Make the plot
  p_516 <- ggplot(eval_516, aes(x=qual_depth, y=percentqPCR, fill=Uses.516.for.Arc)) + 
    geom_boxplot() +
    scale_y_continuous(limits=c(0, 1), labels=percent) +
    scale_fill_manual(name="Uses 516\nfor Archaea", values=c("gray50", "white")) +
    #geom_text(data=p_labels, aes(x=qual_depth+0.5*as.logical(Uses.516.for.Arc), y=0.9, colour=Uses.516.for.Arc, label=n.points)) +
    xlab("depth bin") +
    ylab("percent archaea by qPCR") +
    coord_flip() +
    theme(legend.position="top")
  ggsave("~/Dropbox/Metadata Analysis/Revision for AEM/plots/evaluation_of_516.tiff", height=4, width=clm, units="in", dpi=900, compression="lzw", type="cairo")
  
  list(p_516=p_516,
       med_and_IQR=med_and_IQR,
       t_test_516=t_test_516,
       n_points=n_points)
}