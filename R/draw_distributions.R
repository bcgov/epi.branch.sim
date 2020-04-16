# draw_distributions.R
# Contains the functions that define how various parameters are drawn.
# These functions should be called whenever a stochastic parameter is needed
# so that it is easy swap between different methods. Generally, the
# description for the distributions are defined in `sim.params`

#' Draw incubation periods for new cases from distribution
#' defined in \code{sim.params$incub_params}
#'
#' @param n           Number of cases required
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim.params$incub_params} Current distributions possible are:
#'                    \itemize{
#'                    \item \code{lognormal}: Delay is drawn from lognormal distribution with
#'                    attributes "meanlog" and "sdlog" given in \code{sim.params$incub_params}.
#'                    for traced cases. Untraced cases have a further delay of "delay" days.
#'                    \item \code{weibull}: Delay is drawn from a Weibull distribution with
#'                    attributes "shape" and "scale" given in \code{sim.params$incub_params}.
#'                    Traced cases have their delay set to zero days.
#'                    }
#' @return A vector of length n for case incubation period (double)
draw_incubation_period <- function(n, sim.params){
  if (sim.params$incub_params$dist=='lognormal'){
    return(rlnorm(n, meanlog=sim.params$incub_params$meanlog, sdlog=sim.params$incub_params$sdlog))
  } else if (sim.params$incub_params$dist=='weibull'){
    return(rweibull(n, shape=sim.params$incub_params$shape, scale=sim.params$incub_params$scale))
  }
}

#' Determine whether a new case is symptomatic from a uniform distribution
#'
#' @param n           Number of cases required
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, here,
#'                    the \code{sim.params$p.sym} value is used as the probability of a case
#'                    being symptomatic.
#' @return A boolean vector of length n for whether a case is symptomatic
draw_symptomatic_status <- function(n, sim.params){
 return(runif(n) < sim.params$p.sym)
}

#' Draw onset-to-isolation periods for new cases from distribution
#' defined in \code{sim.params$iso_delay_params}
#'
#' @param state_df    \code{state_df} object for the simulation
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim.params$iso_delay_params}. Current distributions possible are:
#'                    \itemize{
#'                    \item \code{uniform_delay}: Delay is drawn from uniform distribution with
#'                    attributes "min" and "max" given in \code{sim.params$iso_delay_params}.
#'                    for traced cases. Untraced cases have a further delay of "delay" days.
#'                    \item \code{Hellewell}: Delay is drawn from a Weibull distribution with
#'                    attributes "shape" and "scale" given in \code{sim.params$iso_delay_params}.
#'                    Traced cases have their delay set to zero days.
#'                    }
#' @return A vector of length n for case onset-to-isolation period (double)
draw_isolation_delay_period <- function(state_df, sim.params){
  iso_delay_params <- sim.params$iso_delay_params
  if (iso_delay_params$dist=='uniform_delay'){
    n<-nrow(state_df)
    # currently, 2-3 days if traced, 4-5 days if not
    # so, uniformly sample between 2 and 3, then add 2 days if not traced
    iso_delay <- runif(n,min=iso_delay_params$min,max=iso_delay_params$max)
    iso_delay <- iso_delay + (1-state_df$is_traced)*iso_delay_params$delay
  } else if (iso_delay_params$dist=='Hellewell'){
    n<-nrow(state_df)
    # Draw from Weibull distribution initially
    iso_delay<-rweibull(
      n,
      shape=iso_delay_params$shape,
      scale=iso_delay_params$scale
    )
    # Set delay to zero for traced cases
    iso_delay <- iso_delay * (1-state_df$is.traced)
  }
  return(iso_delay)
}

#' Determine whether a case is being traced by public health
#'
#' Draws a uniform random number and compares to trace probability (which can be variable or
#' fixed) to determine whether a case is being traced or not.
#'
#' The probability of a case being traced, \code{p.trace}, can either be a value that varies
#' with cluster size (if \code{sim.params$vary.trace} is \code{TRUE}) or a constant value equal
#' to \code{sim.params$p.trace} if \code{sim.params$vary.trace} is \code{FALSE}.
#'
#' When variable tracing based on cluster size is used, \code{n} is assumed to be the cluster size
#' because it is assumed that this function would be called be \code{create_state_df} for all
#' secondary infections together. The breakpoints in cluster size groupings are defined by
#' \code{c(0,sim.params$cluster_breaks,Inf)} with intervals being closed on the left and open
#' on the right. For instance, if \code{sim.params$cluster_breaks=c(10,25)} then \code{n=9} would
#' be in the first interval, \code{n=10} would be in the second interval and \code{n=25} would be
#' in the third interval.
#'
#' @param n           Number of cases to consider
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, where
#'                    all of the information needed to determine the probability of a case being
#'                    traced (\code{p.trace}) are found within:
#'                    \itemize{
#'                    \item \code{sim.params$vary.trace}: Boolean to determine whether
#'                    \code{p.trace} varies based on cluster size or whether a single constant
#'                    value for \code{p.trace} is used.
#'                    \item \code{sim.params$cluster_breaks}: Defines the breakpoints for
#'                    different \code{p.trace} based on cluster size. Zero and Infinity as the
#'                    lower and upper bounds are implied. Ignored if \code{sim.params$vary.trace}
#'                    is \code{FALSE}.
#'                    \item \code{sim.params$cluster_p.trace}: Defines the \code{p.trace} value for
#'                    each cluster size group defined by \code{sim.params$cluster_breaks}. Length
#'                    of this should be 1 larger than the length of \code{sim.params$cluster_breaks}.
#'                    Ignored if \code{sim.params$vary.trace} is \code{FALSE}.
#'                    \item \code{sim.params$p.trace}: If \code{sim.params$vary.trace} is \code{FALSE}
#'                    then this is the constant \code{p.trace} value used at all times. Ignored if
#'                    \code{sim.params$vary.trace} is \code{TRUE}.
#'                    }
#' @return A vector of length n for case onset-to-isolation period (double)
draw_traced_status <- function(n, sim.params){
  if (sim.params$vary.trace){ # use a variable tracing model
    # Determine which cluster group n is in
    cluster_size_ind <- cut(n,
                            c(0,sim.params$cluster_breaks,Inf),
                            right=FALSE,
                            labels=1:length(sim.params$cluster_p.trace))
    # Assign p_trace based on cluster group size
    p_trace<-sim.params$cluster_p.trace[cluster_size_ind]
    #
    return(runif(n) < p_trace)
  } else{ # otherwise, use fixed tracing value
    return(runif(n) < sim.params$p.trace)
  }

}

#' Determine infection length of new cases
#'
#' Currently a constant number for all cases
#'
#' @param n           Number of cases required
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, here,
#'                    the \code{sim.params$infect.dur} value is used the infection length in days.
#' @return A vector of length n for infection duration (double)
draw_infection_length <- function(n, sim.params){
  return(rep(sim.params$infect.dur,n))
}

#' Assign social distancing behaviour to new cases
#'
#' Input parameters determine the social distancing behaviour for different population groups
#' as well as the probability of a new case being in a given group. Distancing behaviour is
#' defined as a number between 0 and 1, representing the relative number of contacts the
#' subject encounters after distancing begins (i.e. 0.6 means the subject reduces contact
#' to 60% of the normal amount.) There is also a parameter to determine the simulation day
#' number where this behaviour begins.
#'
#' @param n_cases     Number of cases required
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters.
#'                    Here, the sim.params$social_dist_params object contains the information
#'                    needed. This list should have the following entries:
#'                    \itemize{
#'                    \item \code{sd_factors}: A list of social distancing factors with length
#'                    equal to the number of different social distancing groups.
#'                    \item \code{p.group}: A list of probabilities of being in one of the social
#'                    distancing groups. Length should be the same as \code{sd_factors}.
#'                    \item \code{delay}: The simulation day number where social distancing
#'                    begins. All new cases prior to this time will have no distancing factor
#'                    applied (i.e. factor is set to 1).
#'                    }
#' @param sim.status  \code{sim.status} object (a list) containing simulation state vector
#' @param initialize  A boolean indicating whether these cases are the first cases of the
#'                    simulation (sets sd_factor to 1). Defaults to \code{FALSE}.
#' @param import      A boolean indicating whether these new cases are imported cases instead
#'                    of secondary infections. Imported cases do not practice social distancing
#'                    (sd_factor set to 1). Defaults to \code{FALSE}.
#' @return A vector of length n for the sd_factor (double)
draw_sd_factor <- function(n_cases, sim.params, sim.status, initialize=FALSE, import=FALSE){
  sd_params <- sim.params$social_dist_params
  n_groups <- length(sd_params$sd_factors)
  if (initialize){ # no social distancing for initial cases
    sd_factor<-rep(1,n_cases)
  }
  else if (import){ # no social distancing for imported cases
    sd_factor<-rep(1,n_cases)
  }
  else if (sim.status$t < sd_params$delay){ # no social distancing until delay time
    sd_factor<-rep(1,n_cases)
  } # After delay time...
  else if (n_groups==1){
    sd_factor<-rep(sd_params$sd_factors, n_cases)
  }
  else{
    sd_factor <- sample(x=sd_params$sd_factors, size=n_cases, prob=sd_params$p.group, replace=TRUE)
  }
  return(sd_factor)
}

#' Stage secondary infections when new cases are generated
#'
#' When a new case is generated, all the secondary infections are also generated and
#' loaded into the \code{state_df} row for the primary case. The infections are only
#' actually loaded into the simulation at the correct time, through the function
#' \code{generate_secondary_infections()} in generate_new_infections.R. To generate
#' secondary infections, first the number of potential secondary infections is drawn.
#' Then, serial intervals are generated for each potential secondary infection. Finally,
#' potential infections are rejected by user-defined criteria, such as being after the
#' isolation time for the primary case. Currently, only one type of secondary infection
#' algorithm is implemented (the one from Hellewell et al.). The serial intervals are
#' generated by a call to \code{draw_serial_interval}.
#'
#' @param state_df    \code{state_df} object for the newly generated cases
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters.
#'                    Specifically, the \code{sim.params$sec_infect_params} object contains
#'                    the parameters for generating secondary infections.
#'                    This list should have the following entries:
#'                    \itemize{
#'                    \item \code{type}: Defines the algorithm to be used. Currently, only
#'                    "Hellewell" is implemented. This uses a negative binomial distribution
#'                    to draw the number of potential secondary infections with a mean equal
#'                    to R0*sd_factor and a dispersion equal to the value set below. Serial
#'                    intervals are drawn but only potential infections with serial intervals
#'                    earlier than the primary case's time of isolation are kept.
#'                    \item \code{disp}: The dispersion value for the negative binomial distribution
#'                    used to determine the number of potential secondary infections.
#'                    }
#' @param sim.status  \code{sim.status} object (a list) containing simulation state vector
#' @param initialize  A boolean indicating whether these cases are the first cases of the
#'                    simulation, which is passed to \code{draw_sd_factor}.
#'                    Defaults to \code{FALSE}.
#' @param import      A boolean indicating whether these new cases are imported cases instead
#'                    of secondary infections, which is passed to \code{draw_sd_factor}.
#'                    Defaults to \code{FALSE}.
#' @return A list containing
#' \itemize{
#'   \item \code{sd_factor} - A vector of sd_factor for each case (double)
#'   \item \code{n} - A vector with number of accepted secondary infections for each case (int)
#'   \item \code{serial} - A list where each entry corresponds to each case and contains a
#'   vector of serial intervals of accepted secondary infections.
#'   \item \code{non_infects} - A list where each entry corresponds to each case and contains a
#'   vector of serial intervals of rejected secondary infections, for debugging.
#' }
#'
draw_sec_infects_df <- function(state_df, sim.params, sim.status, initialize=FALSE, import=FALSE){
  sec_infect_params <- sim.params$sec_infect_params
  n_cases = nrow(state_df)
  if (sec_infect_params$type=='Hellewell'){
    # Following Hellewell et al, for each case, determine
    # number of sec. infections and serial interval of each sec. infection
    col_names <- c('n_infect', 'serial_int')
    ###sec_infects_df <- data.frame(matrix(nrow=n_cases,ncol=n_cols, dimnames=list(NULL,col_names)))
    # Determine which social distancing group cases belong to
    sd_factor <- draw_sd_factor(n_cases, sim.params, sim.status, initialize=initialize, import=import)
    # Determine number of secondary infections drawn from neg. binomial
    mean_infect <- sim.params$R0 * sd_factor
    disp_infect <- sec_infect_params$disp
    #cat(sprintf("check n_cases = %f,\tmu = %f,\t size=%f \n",n_cases,mean_infect,disp_infect))
    n.infect <- rnbinom(n_cases, mu=mean_infect, size=disp_infect)
    # Determine serial interval of each secondary infection for each infector source case
    if (sim.params$serial_int_params$dist=='gamma'){
      serial.int <- mapply(draw_serial_interval, # FUN to be called on each X
                           n.infect, # X to loop over
                           MoreArgs=list(sim.params=sim.params), # additional required input for FUN
                           SIMPLIFY = FALSE) # force return as list
    } else if (sim.params$serial_int_params$dist=='skewed_norm'){
      serial.int <- mapply(draw_serial_interval, # FUN to be called on each X1, X2
                           n.infect,             # X1 to loop over
                           state_df$incubation_length, #X2 to loop over
                           MoreArgs = list(sim.params = sim.params),# additional required input for FUN
                           SIMPLIFY = FALSE) # force return as list
    }
    # Sort serial intervals in ascending order
    serial.int<-lapply(serial.int,sort)
    ## Determine the last day contagious for each case
    # Switches for each contagious scenario
    is_T_and_S <- state_df$is_traced * state_df$is_symptomatic # traced and symptomatic
    is_T_and_nS <- state_df$is_traced * (1-state_df$is_symptomatic) # traced and not sympt
    is_nT_and_S <- (1-state_df$is_traced) * state_df$is_symptomatic # not traced and sympt
    is_nT_and_nS <- (1-state_df$is_traced) * (1-state_df$is_symptomatic) # not traced and not sympt
    # Length of contagious period for each scenario
    T_and_S_time <- state_df$incubation_length # isolated right after symptoms appear
    T_and_nS_time <- state_df$infection_length # no isolation if no symptoms
    nT_and_S_time <- state_df$incubation_length + state_df$isolation_delay # isolated some time after symptoms
    nT_and_nS_time <- state_df$infection_length # no isolation if no symptoms
    last_day_contagious <- is_T_and_S   * T_and_S_time   +
                           is_T_and_nS  * T_and_nS_time  +
                           is_nT_and_S  * nT_and_S_time  +
                           is_nT_and_nS * nT_and_nS_time
    # Find index of last infection that will happen
    last_infect_ind <- lapply(seq_along(serial.int), function(ii,serial,last_day){
      sum(serial[[ii]]<last_day[ii]) # N.B. serial is sorted!
    }, serial=serial.int, last_day=last_day_contagious)
    # Update number of infections that will occur
    n.infect<-unlist(last_infect_ind) # turn it back into a vector
    # Split the serial interval list.
    # Keep the valid infections for model, store the rest for record keeping
    serial_keep <- lapply(seq_along(serial.int), function(ii, serial, last_ind){
      if (last_ind[ii]==0){
        return(numeric(0))
      }
      else if (last_ind[ii]==length(serial[[ii]])){
        return(serial[[ii]])
      } else {
        return(serial.int[[ii]][1:last_ind[ii]])
      }
    }, serial=serial.int, last_ind=n.infect)
    serial_reject <- lapply(seq_along(serial.int), function(ii, serial, last_ind){
      if (last_ind[ii]==0 | last_ind[[ii]]==length(serial[[ii]])){
        return(numeric(0))
      }
      else {
        return(serial.int[[ii]][(last_ind[ii]+1):length(serial.int[[ii]])])
      }
    }, serial=serial.int, last_ind=n.infect)
    # Will probably move this outside of the if-block when other methods added
    return(list(sd_factor=sd_factor,n=n.infect,serial=serial_keep,non_infects=serial_reject))
  }
}

#' Generate serial intervals for potential secondary infections
#'
#' Serial intervals for each case's potential secondary infections are drawn from
#' distribution defined in \code{sim.params$serial_int_params}
#'
#' @param n                   Number of potential secondary infections for this case
#' @param incubation_length   Incubation period of this case (used in skewed normal distribution only)
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim.params$serial_int_params}. Current distributions possible are:
#'                    \itemize{
#'                    \item \code{dummy}: Uniform number between 1 and 5 (used for testing only)
#'                    \item \code{skewed_norm}: Drawn from skewed normal distribution with xi set
#'                    to the case's incubation period, and attributes "omega" and "alpha" given
#'                    in \code{sim.params$serial_int_params}. This is what Hellewell uses.
#'                    \item \code{gamma}: Drawn from a gamma distribution with attributes
#'                    "shape" and "rate" given in \code{sim.params$serial_int_params}.
#'                    }
#' @return A vector of length n for serial intervals of the case's potential secondary infections
draw_serial_interval <- function(n, incubation_length, sim.params){
  if (sim.params$iso_delay_params$dist=='dummy'){
    return(runif(n,min=1,max=5))
  } else if (sim.params$serial_int_params$dist=='skewed_norm'){
    sn_xi = incubation_length # case incubation period
    sn_omega = sim.params$serial_int_params$omega
    sn_alpha = sim.params$serial_int_params$alpha
    serial_ints <- sn::rsn(n, xi=sn_xi, omega=sn_omega, alpha=sn_alpha)
    return(as.numeric(serial_ints))
  } else if (sim.params$serial_int_params$dist=='gamma'){
    serial_ints <- rgamma(n, shape=sim.params$serial_int_params$shape, rate=sim.params$serial_int_params$rate)
    return(serial_ints)
  }
}
