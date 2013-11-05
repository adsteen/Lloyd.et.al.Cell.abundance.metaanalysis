##' Reproduces analysis from Lloyd et al (submitted)
##' 
##' @description Reproduces entire analysis from Lloyd et al, 2013, Applied and Environmental Microbiology. 
##' @param print_plots If true, print the plots to screen
##' @param save_plots If true, save plots. There is currently no option to change the filenames of the saved plots.
##' @param reload_data If true, reload the data from .xlsx files. Otherwise read data from .Rdata files. This is passed along to `read_data`. 
##' @param fast_calc If true, skip the bootstrapped parts of the analysis, which take the lion's share of computational time.
##' @param n_reps Number of reps for bootstrapped analyses.
##' @return Returns the three central data frames used in the analysis: `all_data`, `corrected_seds`, and `corrected_sw`, as named elements in a list. Note that `corrected_seds` and `corrected_sw` are, collectively, redundant with `all_data`
##' @export
reproduce_research <- function(print_plots=TRUE, save_plots=FALSE, reload_data=FALSE, fast_calc=TRUE, n_reps=1000) {
  
  ###
  # Uncomment these when running this function as a script
  ###
#   reload_data <- FALSE
#   print_plots <- TRUE
#   save_plots <- FALSE
#   fast_calc <- TRUE
#   n_reps <- 100
#   point_size <- 0.5
  
  #====================
  # global figure settings
  #====================
  
  clm <- 3.42 # Single column width (inches)
  clm2 <- 7.08 # Double column width (inches)
  myDPI <- 900 # Resolution (DPI)
  
  theme_set(theme_bw() + 
              theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_rect(colour="black"),
                    text=element_text(size=8)))
  
  #======================
  # Load data: 
  #     -from .Rdata file supplied with package (by default), 
  #     -from other .Rdata files (by supplying filenames) 
  #     -or from .xlsx files (reload.from.xlsx=TRUE)
  #    
  #=====================
  
  data_list <- read_data(reload=reload_data)
  all_data <- data_list[["all_data"]] 
  corrected_seds <- data_list[["corrected_seds"]]
  corrected_sw <- data_list[["corrected_sw"]]
  
  
  ##############
  # Figure 1: total prokaryotes, by CARDFISH, vs. total cells, by direct counts, in sw and sediments
  ##############
  p_cells_vs_fish <- plot_cell_vs_fish(all_data)
  
  if (print_plots) {
    print(p_cells_vs_fish)
  }
  if (save_plots) {
    tiff("plots/p_cells_vs_fish.tif", height=3, width=clm, units="in", res=myDPI, compression="lzw", type="cairo")
    print(p_cells_vs_fish)
    dev.off()
  }
  
  
  ################
  # Fig 2: FISH / CARDFISH yield in sediments, by permeabilization method and core
  ###############
  
  # Parameters needed by the figure
  yaxsize <- 7
  yield_label <- expression(frac("Bacteria + Archaea (FISH or CARD-FISH)", "Total Cell Count"))
  colors <- brewer.pal(n=5, name="RdBu")
  
  # Calculate the median and interquartile range of seawater FISH yields
  sw_quant <- quantile(corrected_sw$FISH.yield, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)
  
  # Boxplot of yields for each core, colored by permeabilization method
  p_yield_by_core <- yield_by_core(corrected_seds=corrected_seds, 
                                   sw_quant=sw_quant,
                                   yield_label=yield_label)
  
  if (print_plots) {
    print(p_yield_by_core)
  }
  
  if (save_plots) {
    tiff("plots/yield_by_core.boxplot.tif" , height=5, width=clm2, units="in", res=900, compression="lzw", type="cairo")
    print(p_yield_by_core)
    dev.off()
  }
  
  all_perm_boxplots <- boxplots_by_perm(corrected_seds, sw_quant)
  
  vp1 <- viewport(x=0.25, y=0.75, height=0.48, width=0.5)
  vp2 <- viewport(x=0.75, y=0.75, height=0.48, width=0.5)
  vp3 <- viewport(x=0.25, y=0.25, height=0.48, width=0.5)
  vp4 <- viewport(x=0.75, y=0.25, height=0.48, width=0.5)
  
  if (print_plots) {
    print(all_perm_boxplots$p_yield_all, vp=vp1)
    print(all_perm_boxplots$p_yield_shallow, vp=vp2)
    print(all_perm_boxplots$p_percent_arc_all_depths, vp=vp3)
    print(all_perm_boxplots$p_percent_arc_shallow, vp=vp4)
  }
  
  if (save_plots) {
    #print the actual plots
    tiff("plots/boxplots_percent_arc_by_perm_method.tif", height=5, width=clm2, units="in", res=myDPI, compression="lzw", type="cairo")
    print(all_perm_boxplots$p_yield_all, vp=vp1)
    print(all_perm_boxplots$p_yield_shallow, vp=vp2)
    print(all_perm_boxplots$p_percent_arc_all_depths, vp=vp3)
    print(all_perm_boxplots$p_percent_arc_shallow, vp=vp4)
    dev.off()
  }
  
  
  ######
  # Fig 4: qPCR quantifications
  ######
  # This will throw multiple warnings; that's OK
  p_qPCR_list <- make_qPCR_plots(all_data, corrected_seds)
  
  vp_left <- viewport(x=0.25, y=0.5, height=1, width=0.5)
  vp_top_right <- viewport(x=0.75, y=0.75, height=0.5, width=0.5)
  vp_bottom_right <- viewport(x=0.75, y=0.25, height=0.5, width=0.5)
  
  if(print_plots) {
    print(p_qPCR_list$p_qPCRtotal_v_totalcells, vp=vp_left)
    print(p_qPCR_list$p_qPCR_sum_bac_arc_v_universal, vp=vp_bottom_right)
    print(p_qPCR_list$p_var_v_var, vp=vp_top_right)
  }
  
  if (save_plots) {
    tiff("plots/all_qPCR_plots_no_legend.tif", height=7, width=clm2, 
         units="in", res=900, compression="lzw", type="cairo")
    print(p_qPCR_list$p_qPCRtotal_v_totalcells, vp=vp_left)
    print(p_qPCR_list$p_qPCR_sum_bac_arc_v_universal, vp=vp_bottom_right)
    print(p_qPCR_list$p_var_v_var, vp=vp_top_right)
    dev.off()
  }
  
  # Print fig 3a, again, with a legend, for supplemental figures
  p_qPCR_supplemental <- p_qPCR_list$p_qPCRtotal_v_totalcells +
    theme(legend.position="bottom")
  if (print_plots) {
    print(p_qPCR_supplemental)
  }
  
  if (save_plots) {
    ggsave("plots/qPCR_plot_for_supplemental.png", p_qPCR_supplemental, height=8, width=clm*1.5, 
           units="in", dpi=myDPI)
  }
  
  
  ##########
  # Fig 5: Trends in archaea & bacteria (relative & absolute) by depth 
  ##########
  
  # Seawater depth trends
  sw_profile_list <- sw_depth_profiles(corrected_sw)
  
  # Sediment depth trends
  sed_percent_arc_v_depth_plot <- sed_percent_arc_v_depth(corrected_seds, point_size=point_size)
  sed_bac_arc_v_depth_plots <- sed_bac_and_arc_v_depth(corrected_seds, point_size=point_size)
  
  h_top <- 0.44
  h_bottom <- 1-h_top
  y_top <- h_bottom + h_top/2
  y_bottom <- h_bottom/2
  vp1 <- viewport(height=h_top, width=1/3, x=1/6, y=y_top)
  vp2 <- viewport(height=h_top, width=1/3, x=3/6, y=y_top)
  vp3 <- viewport(height=h_top, width=1/3, x=5/6, y=y_top)
  vp4 <- viewport(height=h_bottom, width=1/3, x=1/6, y=y_bottom)
  vp5 <- viewport(height=h_bottom, width=1/3, x=3/6, y=y_bottom)
  vp6 <- viewport(height=h_bottom, width=1/3, x=5/6, y=y_bottom)
  
  if (print_plots) {
    print(sw_profile_list$p_bac_profile, vp=vp1)
    print(sw_profile_list$p_arc_profile, vp=vp2)
    print(sw_profile_list$p_percent_profile, vp=vp3)
    print(sed_bac_arc_v_depth_plots$p_sed_bac, vp=vp4)
    print(sed_bac_arc_v_depth_plots$p_sed_arc, vp=vp5)
    print(sed_percent_arc_v_depth_plot$p_sed_frac_arc_v_depth, vp=vp6)
  }
  
  if (save_plots) {
    tiff("plots/all_depth_trends_segmented.tif", height=6, width=clm2, 
         units="in", res=900, compression="lzw", type="cairo")
    print(sw_profile_list$p_bac_profile + theme(plot.margin=unit(c(1, 1, 1, 1), "mm")), vp=vp1)
    print(sw_profile_list$p_arc_profile + theme(plot.margin=unit(c(1, 1, 1, 1), "mm")), vp=vp2)
    print(sw_profile_list$p_percent_profile + theme(plot.margin=unit(c(1, 2, 1, 1), "mm")), vp=vp3)
    print(sed_bac_arc_v_depth_plots$p_sed_bac + theme(plot.margin=unit(c(1, 1, 1, 1), "mm")), vp=vp4)
    print(sed_bac_arc_v_depth_plots$p_sed_arc + theme(plot.margin=unit(c(1, 1, 1, 1), "mm")), vp=vp5)
    print(sed_percent_arc_v_depth_plot$p_sed_frac_arc_v_depth + 
            theme(plot.margin=unit(c(1, 2, 5, 1), "mm"),
                    axis.title.x=element_text(vjust=-1)), vp=vp6)
    dev.off()
  }
  
  
  ##############
  # Miscellaneous analysis that goes into the text 
  #    Note: this should be run after all the figure functions, as it relies on some of the 
  #          data frames generated by those functions
  ##############
  
  ################
  # Comparisons of yield by different methodology
  ################
  
  if (fast_calc == FALSE) {
    system.time(sw_yield_plot_list <- make_sw_yield_boxplots(corrected_sw, n=n_reps))
    
    if (print_plots) {
      # Print the plots
      grid.arrange(sw_yield_plot_list[[1]], 
                   sw_yield_plot_list[[2]],
                   sw_yield_plot_list[[3]],
                   sw_yield_plot_list[[4]],
                   sw_yield_plot_list[[5]],
                   main="seawater")
      
      # Print the summary statistics
      sw_yield_summary_statistics <- llply(sw_yield_plot_list[6:10], print)
      save(sw_yield_summary_statistics, file="data/sw_yield_summary_statistics.RData")
    }
    
    if (save_plots) {
      # Save the plots
      tiff("plots/sw_yield_boxplots_permutation_stats.tif", width=clm2, height=clm2, units="in",res=myDPI, compression="lzw")
      grid.arrange(sw_yield_plot_list[[1]], 
                   sw_yield_plot_list[[2]],
                   sw_yield_plot_list[[3]],
                   sw_yield_plot_list[[4]],
                   sw_yield_plot_list[[5]],
                   main="seawater")
      dev.off()
      
      # Print the summary statistics to .pdf tables
      fn_vec_sw <- paste("tables/", 
                       c("sw_hyb_method", "sw_stain", "sw_automatic", "sw_arc_probe", "sw_bac_probe"), 
                       ".pdf", sep="")
      for (i in 6:10) {
        pdf(fn_vec_sw[i-5])
        grid.table(sw_yield_plot_list[[i]])
        dev.off()
      }
    }
    
    
    # Make sediment plots
    system.time(sed_yield_plot_list <- make_sed_yield_boxplots(corrected_seds, n=n_reps))
    
    if (print_plots) {
      grid.arrange(sed_yield_plot_list[[1]],
                   sed_yield_plot_list[[2]],
                   sed_yield_plot_list[[3]],
                   main="sediments")
      l_ply(sed_yield_plot_list[4:6], print)
    }
    
    if (save_plots) {
      tiff("plots/sed_yield_boxplots.tif", height=8, width=clm, units="in", res=myDPI, compression="lzw", type="cairo")
      grid.arrange(sed_yield_plot_list[[1]],
                   sed_yield_plot_list[[2]],
                   sed_yield_plot_list[[3]],
                   main="sediments")
      dev.off()
      
      # Save tables of summary statistics
      fn_vec_sed <- paste("tables/", 
                          c("sed_stain", "sed_bac_probe", "sed_hyb_method"), 
                          ".pdf", sep="")
      for (i in 4:6) {
        pdf(fn_vec_sed[i-3])
        grid.table(sed_yield_plot_list[[i]])
        dev.off()
      }
    }
  }
  
  
  
  ############
  # Yields for intertidal sediments
  ############
  p_intertidal_yield <- intertidal_yield_fig(corrected_seds, corrected_sw)
  
  if (print_plots) {
    print(p_intertidal_yield)
  }
  
  if (save_plots) {
    tiff("plots/supp_fig_intertidal_sediments_yield.tiff", height=4, width=clm2, units="in", res=myDPI, type="cairo")
    print(p_intertidal_yield)
    dev.off()
  }
  
  #########
  # Evaluation of 516 as qPCR primer
  #########
  qPCR_516 <- qPCR_516_evaluation(corrected_seds)
  
  if(print_plots) {
    print(qPCR_516$p_516)
    print(qPCR_516$med_and_IQR)
  }
  
  if (save_plots) {
    ggsave("plots/evaluation_of_516.tiff", height=4, width=clm,units="in", res=myDPI, type="cairo")
    pdf("tables/evaluation_of_516.pdf")
    grid.table(qPCR_516$med_and_IQR)
    dev.off()
  }
  
  
  #############
  # Fig S3: Absolute quantities of bacteria & Archaea, in sw and seds, quantified by best practices
  ############

  #=========
  # Return all data for later use
  #=========
  
  return_list <- list(all_data=all_data,
                      corrected_sw=corrected_sw,
                      corrected_seds=corrected_seds)
  return_list
}