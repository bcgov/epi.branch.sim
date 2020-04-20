# plotting_functions.R
# Function to make plots!

library(ggplot2)
library(ggpubr)
library(RColorBrewer) # also needs scales

#' Generates cumulative and daily counts plot object
#'
#' This function creates a two-panel figure. The top panel shows the cumulative case
#' counts and the bottom shows the daily new case counts. The x-axis for both is the
#' date. The required input parameter, \code{results} is a \code{tibble} with the
#' following columns:
#' \itemize{
#'   \item date (\code{\link[base]{Dates}} objects),
#'   \item n.iso (integer, total number of isolated cases on given day)
#'   \item dn.iso (integer, number of new isolated cases on a given day)
#'   \item reported (integer, number of total reported cases on a given day)
#'   \item reported_daily (integer, number of new reported cases on a given day)
#' }
#'
#' The dataframe that is returned by \code{\link{run_scenarios}} can be
#' used as a starting point to build the \code{results} tibble needed for input to this
#' function. The \code{day} column from the \code{\link{run_scenarios}} output is the
#' day number in the simulation and can be converted to a \code{\link[base]{Dates}} object
#' with a real calendar date using
#' \preformatted{
#'  results <- as_tibble(results) \%>\%
#'    mutate(date=lubridate::as_date("2020-03-01")+lubridate::days(day)) \%>\%
#'    mutate_at(c("R0","initial_n","p.trace","p.symp"),as.factor)
#' }
#' where March 1, 2020 is simulation day 0 in this example. The column \code{dn.iso} can be
#' computed from \code{n.iso} using
#' \preformatted{
#' dn.iso <- diff(results$n.iso) # the first entry (day 0) will be negative
#' dn.iso <- c(0,dn.iso) # diff returns length n-1, add 0 to first position
#' dn.iso[dn.iso<0] <- 0 # set the day 0 to zero
#' results$dn.iso <- dn.iso
#' }
#' Finally, if using \code{plot_reported=TRUE}, then to add the \code{reported} column,
#' it would be a good idea to use \code{NA} for the dates you don't want to plot as well
#' as future dates where reported counts are not available. The \code{reported_daily}
#' column can be computed in a similar way to the \code{dn.iso} example above.
#'
#' @param results         A tibble with required columns described here
#' @param start_date      Start of x-axis, string in "YYYY-MM-DD" format
#' @param end_date        End of x-axis, string in "YYYY-MM-DD" format
#' @param title_string    Title of plot
#' @param max_y_cum,max_y_daily   Upper limit on y-axis of cumulative and daily plots,
#'                                defaults to \code{NULL}, which sets limit based on data
#' @param log_y_cum,log_y_daily   If \code{TRUE}, uses a log scale for y-axis on
#'                                cumulative and daily plots. Defaults to \code{FALSE}
#' @param plot_reported   Plots reported counts if \code{TRUE}. If used, \code{results}
#'                        must have \code{reported} and \code{reported_daily} as columns.
#' @param label_projected Label for curves showing model predicted case counts
#' @param label_reported  Label for curves showing reported case counts
#' @param dark_quantiles  Quantile limits to use for darker shade. Default: c(0.25,0.75)
#' @param light_quantiles Quantile limits to use for lighter shade. Default: c(0.05,0.95)
#' @param legend_pos      Vector for legend position, relative to plot axis. Default: c(0.25,0.75)
#' @return A ggplot object with two individual plots showing cumulative case counts (top) and
#'         daily new case counts (bottom) as a function of time.
plot_iso_case_counts <- function(results, start_date, end_date,
                                 title_string='',
                                 max_y_cum=NULL, max_y_daily=NULL,
                                 log_y_cum=FALSE, log_y_daily=FALSE,
                                 plot_reported=FALSE,
                                 label_projected='Projected counts',
                                 label_reported='Reported counts',
                                 dark_quantiles=c(0.25,0.75),
                                 light_quantiles=c(0.05,0.95),
                                 legend_pos=c(0.25,0.75)){
  paired.cols <- brewer.pal(12, "Paired")
  colours <- c(label_projected = paired.cols[2], label_reported = paired.cols[8])

  # Create cumulative counts plot object with median and quantile bands
  p1 <- ggplot(results, aes(x=date)) +
    stat_summary(geom="ribbon", aes(y=n.iso),
                 fun.ymin = function(x) quantile(x, light_quantiles[1]),
                 fun.ymax = function(x) quantile(x, light_quantiles[2]),
                 fill=paired.cols[1], alpha=0.8) +
    stat_summary(geom="ribbon", aes(y=n.iso),
                 fun.ymin = function(x) quantile(x, dark_quantiles[1]),
                 fun.ymax = function(x) quantile(x, dark_quantiles[2]),
                 fill=paired.cols[2], alpha=0.5) +
    stat_summary(geom="line", aes(y=n.iso, color=label_projected),
                 fun.y=median, size=1.2)
  # Add reported counts if set
  if (plot_reported){
    p1 <- p1 + geom_point(data=results,aes(y=reported, color=label_reported), alpha=0.5)
  }
  # Set plot scales and labels
  p1 <- p1 + xlim(as.Date(start_date),as.Date(end_date)) +
    labs(x="Date", y = "Cumulative cases", color="Cases") +
    ggtitle(title_string) +
    scale_color_manual(values = colours) +
    scale_y_continuous(labels = scales::comma) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          legend.position = legend_pos,
          legend.title = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background = element_rect(colour ="transparent", fill='transparent'),
          axis.title.x=element_blank())
  # Use a log y scale if set
  if (log_y_cum){
    p1 <- p1 + scale_y_log10(labels=scales::comma)
  }

  # Use a user-defined max y value (for comparison across models)
  if (!is.null(max_y_cum)){
    p1 <- p1 + coord_cartesian(ylim = c(0, max_y_cum))
  }

  # Create daily counts plot object with median and quantile bands
  p2 <- ggplot(results,aes(x=date)) +
    stat_summary(geom="ribbon", aes(y=dn.iso),
                 fun.ymin = function(x) quantile(x, light_quantiles[1]),
                 fun.ymax = function(x) quantile(x, light_quantiles[2]),
                 fill=paired.cols[1], alpha=0.8) +
    stat_summary(geom="ribbon", aes(y=dn.iso),
                 fun.ymin = function(x) quantile(x, dark_quantiles[1]),
                 fun.ymax = function(x) quantile(x, dark_quantiles[2]),
                 fill=paired.cols[2], alpha=0.5) +
    stat_summary(geom="line", aes(y=dn.iso),
                 fun.y=median, color=paired.cols[2], size=1.2)
  # Add reported counts if set
  if (plot_reported){
    p2 <- p2 + geom_point(data=results,aes(y=reported_daily, color=label_reported), alpha=0.5)
  }
  # Add plot scales and labels
  p2 <- p2 + xlim(as.Date(start_date),as.Date(end_date)) +
    labs(x="Date", y = "Daily new cases", color="Cases") +
    scale_colour_manual(values = colours) +
    theme(axis.title.x=element_blank(),
          legend.position = "none")
  # Use a log y scale if set
  if (log_y_daily){
    p2 <- p2 + scale_y_log10(labels=scales::comma)
  }
  # Use a user-defined max y value (for comparison across models)
  if (!is.null(max_y_daily)){
    p2 <- p2 + coord_cartesian(ylim = c(0, max_y_daily))
  }

  # Stack p1 and p2 vertically
  p<-ggarrange(p1,p2, ncol=1, nrow=2,
               heights = c(1.2, 0.7), align='v')
  return(p)
}

#' Generates plots of new secondary infections and effective R
#'
#' This function creates a two-panel figure. The top panel shows the daily number of
#' new secondary infections and the bottom shows the effective R. The x-axis for both is the
#' date. The required input parameter, \code{results} is a \code{tibble} with the
#' following columns:
#' \itemize{
#'   \item date (\code{\link[base]{Dates}} objects),
#'   \item n.new_S (integer, total number of secondary cases on given day)
#'   \item Reff (double, estimate of effective R up to a given day)
#' }
#'
#' The dataframe that is returned by \code{\link{run_scenarios}} can be
#' used as a starting point to build the \code{results} tibble needed for input to this
#' function. The \code{day} column from the \code{\link{run_scenarios}} output is the
#' day number in the simulation and can be converted to a \code{\link[base]{Dates}} object
#' with a real calendar date using
#' \preformatted{
#'  results <- as_tibble(results) \%>\%
#'    mutate(date=lubridate::as_date("2020-03-01")+lubridate::days(day)) \%>\%
#'    mutate_at(c("R0","initial_n","p.trace","p.symp"),as.factor)
#' }
#' where March 1, 2020 is simulation day 0 in this example.
#'
#' @param results         A tibble with required columns described here
#' @param start_date      Start of x-axis, string in "YYYY-MM-DD" format
#' @param end_date        End of x-axis, string in "YYYY-MM-DD" format
#' @param title_string    Title of plot
#' @param max_y_newS,max_y_Reff   Upper limit on y-axis of new secondary and Reff plots,
#'                                defaults to \code{NULL}, which sets limit based on data
#' @param log_y_newS      If \code{TRUE}, uses a log scale for y-axis on upper panel (number
#'                        of new secondary cases). Defaults to \code{FALSE}
#' @param dark_quantiles  Quantile limits to use for darker shade. Default: c(0.25,0.75)
#' @param light_quantiles Quantile limits to use for lighter shade. Default: c(0.05,0.95)
#' @return A ggplot object with two individual plots showing cumulative case counts (top) and
#'         daily new case counts (bottom) as a function of time.
plot_Reff <- function(results, start_date, end_date,
                      title_string='',
                      max_y_newS=NULL, max_y_Reff=NULL,
                      log_y_newS=FALSE,
                      dark_quantiles=c(0.25,0.75),
                      light_quantiles=c(0.05,0.95)){
  paired.cols <- brewer.pal(12, "Paired")

  # Generate plot, lines and bands
  p1 <- ggplot(results, aes(x=date)) +
    stat_summary(geom="ribbon", aes(y=n.new_S),
                 fun.ymin = function(x) quantile(x, light_quantiles[1]),
                 fun.ymax = function(x) quantile(x, light_quantiles[2]),
                 fill=paired.cols[7], alpha=0.8) +
    stat_summary(geom="ribbon", aes(y=n.new_S),
                 fun.ymin = function(x) quantile(x, dark_quantiles[1]),
                 fun.ymax = function(x) quantile(x, dark_quantiles[2]),
                 fill=paired.cols[8], alpha=0.5) +
    stat_summary(geom="line", aes(y=n.new_S),
                 color=paired.cols[8], fun.y=median, size=1.2)
  # Set plot scales and labels
  p1 <- p1 + xlim(as.Date(start_date),as.Date(end_date)) +
    labs(x="Date", y = "New secondary infections") +
    ggtitle(title_string) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          axis.title.x=element_blank())
  # Use a log y scale if set
  if (log_y_newS){
    p1 <- p1 + scale_y_log10(labels=scales::comma)
  }
  # Use a user-defined max y value (for comparison across models)
  if (!is.null(max_y_newS)){
    p1 <- p1 + coord_cartesian(ylim = c(0, max_y_newS))
  }

  # Create plot lines and bands
  p2 <- ggplot(results,aes(x=date)) +
    stat_summary(geom="ribbon", aes(y=Reff),
                 fun.ymin = function(x) quantile(x, light_quantiles[1]),
                 fun.ymax = function(x) quantile(x, light_quantiles[2]),
                 fill=paired.cols[1], alpha=0.8) +
    stat_summary(geom="ribbon", aes(y=Reff),
                 fun.ymin = function(x) quantile(x, dark_quantiles[1]),
                 fun.ymax = function(x) quantile(x, dark_quantiles[2]),
                 fill=paired.cols[2], alpha=0.5) +
    stat_summary(geom="line", aes(y=Reff),
                 fun.y=median, color=paired.cols[2], size=1.2)
  # Set plot scales and labels
  p2 <- p2 + xlim(as.Date(start_date),as.Date(end_date)) +
    labs(x="Date", y = "Effective R") +
    theme(axis.title.x=element_blank())
  # Use a user-defined max y value (for comparison across models)
  if (!is.null(max_y_Reff)){
    p2 <- p2 + coord_cartesian(ylim = c(0, max_y_Reff))
  }

  # Stack p1 and p2 vertically
  p<-ggarrange(p1,p2, ncol=1, nrow=2,
               heights = c(1.2, 0.7), align='v')
  return(p)
}
