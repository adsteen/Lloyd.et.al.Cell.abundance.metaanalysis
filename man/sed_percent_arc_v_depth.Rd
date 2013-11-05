\name{sed_percent_arc_v_depth}
\alias{sed_percent_arc_v_depth}
\title{Make depth profiles and calculate linear models for sediment percent Archaea vs depth}
\usage{
  sed_percent_arc_v_depth(corrected_seds,
    point_size = 0.75)
}
\arguments{
  \item{corrected_seds}{Data frame of sediments data}

  \item{point_size}{point size for the saved plots}
}
\value{
  Returns a list with the following elements

  * p_sed_frac_arc_v_depth, a plot of the fraction of
  archaea vs depth in sediments *
  m_sed_frac_arc_depth_qPCR, linear model of fraction of
  archaea vs depth by qPCR * sed_frac_CARD_all_stats,
  correlation coefficients & other statistics from linear
  models of CARD-FISH counts vs depth in sediments *
  sed_frac_est_all_stats, correlation coefficients & other
  statistics from linear models of fraction of Archaea vs
  depth in sediments * frac_est_AIC, data frame of Akaike
  Inforamtion Criterion results from segmented and
  unsegmented linear models * frac_est_break, breakpoint
  (depth, m) in fraction of Archaea by qPCR
}
\description{
  `sed_percent_arc_v_depth`
}

