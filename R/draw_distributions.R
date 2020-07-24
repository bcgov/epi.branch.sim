# Copyright 2020 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.

# draw_distributions.R
# Contains the functions that define how various parameters are drawn.
# These functions should be called whenever a stochastic parameter is needed
# so that it is easy swap between different methods. Generally, the
# description for the distributions are defined in `sim_params`

#' Draw incubation periods for new cases from distribution
#' defined in \code{sim_params$incub_params}
#'
#' @param n           Number of cases required
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim_params$incub_params} Current distributions possible are:
#'                    \itemize{
#'                    \item \code{lognormal}: Delay is drawn from lognormal distribution with
#'                    attributes "meanlog" and "sdlog" given in \code{sim_params$incub_params}.
#'                    for traced cases. Untraced cases have a further delay of "delay" days.
#'                    \item \code{weibull}: Delay is drawn from a Weibull distribution with
#'                    attributes "shape" and "scale" given in \code{sim_params$incub_params}.
#'                    Traced cases have their delay set to zero days.
#'                    }
#' @return A vector of length n for case incubation period (double)
draw_incubation_period <- function(n, sim_params){
  if (sim_params$incub_params$dist=='lognormal'){
    return(rlnorm(n, meanlog=sim_params$incub_params$meanlog, sdlog=sim_params$incub_params$sdlog))
  } else if (sim_params$incub_params$dist=='weibull'){
    return(rweibull(n, shape=sim_params$incub_params$shape, scale=sim_params$incub_params$scale))
  }
}

#' Determine whether a new case is symptomatic from a uniform distribution
#'
#' @param n           Number of cases required
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, here,
#'                    the \code{sim_params$p_sym} value is used as the probability of a case
#'                    being symptomatic.
#' @return A boolean vector of length n for whether a case is symptomatic
draw_symptomatic_status <- function(n, sim_params){
 return(runif(n) < sim_params$p_sym)
}

#' Draw delay to isolation periods for new cases
#'
#' The number of days between symptom onset and isolation, which may be negative if isolation
#' occurs prior to symptom onset (as sometimes the case with traced cases).
#'
#' Within the list object \code{sim_params$iso_delay_params}, the user can define several delay
#' periods based on tracing type, tracing status, and whether the case is practicing distancing.
#' Traced cases have their delays measured from the index case's isolation time, so traced cases
#' may isolate prior to their own symptom onset. Untraced cases have delays measured from the start
#' of their own symptom onset. The untraced timeline is always considered for traced cases, so that
#' if a traced case would have been isolated earlier just because of their symptom onset timeline,
#' they would isolate at this earlier time.
#'
#' Although this function returns the number of days between symptom onset and isolation, the delay
#' returned by this function may be negative if isolation occurs prior to symptom onset.
#'
#' @param state_df    \code{state_df} object for the simulation
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim_params$iso_delay_params}. Current distributions possible are:
#'                    \itemize{
#'                    \item \code{uniform}: Delay is drawn from uniform distributions with
#'                    attributes "min" and "max" given in \code{sim_params$iso_delay_params} for
#'                    each type of case. Min/Max values to be provided are:
#'                       \itemize{
#'                       \item \code{iso_delay_traced_[min|max]}: Range for manually traced cases
#'                       \item \code{iso_delay_untraced_[min|max]}: Range for untraced cases from the
#'                       non-distancing population.
#'                       \item \code{iso_delay_untraced_pd_[min|max]}: Rnage for untraced cases from
#'                       the distancing population (any case with contact_rate less than 1).
#'                       }
#'                    Cases traced by the app (require both index and secondary cases to be app users
#'                    *and* for the secondary case to be app-compliant) have a zero day delay.
#'                    \item \code{Hellewell}: Delay is drawn from a Weibull distribution with
#'                    attributes "shape" and "scale" given in \code{sim_params$iso_delay_params}.
#'                    Traced cases have their delay set to zero days.
#'                    }
#' @param primary_state_df  The \code{state_df} object for the index/primary cases. Defaults to \code{NULL}.
#'                          Only required for traced cases (i.e. not needed when generating initial or
#'                          imported cases).
#' @param primary_case_ids  A list of case_ids for the primary cases. Not required for initial or imported
#'                          cases. Defaults to \code{NULL}.
#' @return A vector of length n for case delay to isolation, measured from start of symptom onset (double)
draw_isolation_delay_period <- function(state_df, sim_params,
                                        primary_state_df=NULL,
                                        primary_case_ids=NULL){
  iso_delay_params <- sim_params$iso_delay_params
  if (iso_delay_params$dist=='uniform'){
    n<-nrow(state_df)
    # Draw an isolation delay time assuming everyone is untraced at first,
    iso_delay <- runif(n,
                       min = iso_delay_params$untraced_min,
                       max = iso_delay_params$untraced_max)
    # Check if primary case info provided (only needed for traced cases)
    if (!is.null(primary_state_df)){
      # Get rows of primary_state_df matching primary cases
      infector_df <- subset(primary_state_df,case_id %in% primary_case_ids)
      # Calculate how many days between now and infector being isolated
      # (may be negative if infector already isolated)
      # HOWEVER, isolated cases shouldn't be causing infections so unlikely to be negative
      infector_df$d_until_iso <- infector_df$incubation_length +
        infector_df$isolation_delay -
        infector_df$days_infected
      # Modify isolation delay based on infector's isolation time, if traced OR
      # if untraced & distancing
      iso_delay2 <- vapply(X=seq_along(iso_delay), FUN=function(ii){
        if (
          state_df[ii,]$is_traced ||
          (state_df[ii,]$is_traced_by_app && state_df[ii,]$is_trace_app_comply)
        ){
          # Get infector isolation time
          d_until_primary_iso <- subset(infector_df,case_id==primary_case_ids[ii])$d_until_iso
          # Get this case's symptom onset time
          d_until_onset <- state_df$incubation_length[ii]
          # Check if new case is traced by app and set appropriate trace delay values
          if (state_df[ii,]$is_traced_by_app && state_df[ii,]$is_trace_app_comply){
            # No delay in tracing, isolate secondary case at infector's isolation
            trace_delay <- 0
          } else {
            # Draw a delay from the manual traced delay distribution
            trace_delay <- runif(1,
                                 min = iso_delay_params$traced_min,
                                 max = iso_delay_params$traced_max)
          }
          # Rest of code/model expects iso_delay to be measured from symptom onset time
          # But we want the traced_delay value to represent the days from the infector's
          # isolation time. So we subtract off this case's symptom onset time and add
          # the number of days until the primary is isolated.
          # This can be negative, meaning we isolate before symptom onset
          # SEE ALSO: draw_sec_infects() in this file
          modified_delay <- trace_delay - d_until_onset + d_until_primary_iso
          # Make sure the modified delay isn't somehow longer than the untraced delay
          # (i.e. if Case A takes a very long time to isolate, Case B isn't traced yet
          # anyhow, so it should just be treated as untraced and isolated at that point)
          return( min(modified_delay, iso_delay[ii]) ) # NB: return is to vapply
        } else{ # If not traced, then only update value if case is also distancing
          if (state_df[ii,]$contact_rate < 1){ # the case is distancing
            return(
              runif(1,min = iso_delay_params$untraced_pd_min,max = iso_delay_params$untraced_pd_max)
            )
          } # not distancing, return original value without modification
          return(iso_delay[ii]) # NB: return is to vapply
        }
      }, FUN.VALUE=999)
      # Update iso_delay vector to account for tracing
      iso_delay <- iso_delay2
    }
  } else if (iso_delay_params$dist=='Hellewell'){
    n<-nrow(state_df)
    # Draw from Weibull distribution initially
    iso_delay<-rweibull(
      n,
      shape=iso_delay_params$shape,
      scale=iso_delay_params$scale
    )

    # Modify actual delay time based on tracing and status of primary case
    # Only when primary_state_df provided (i.e. a secondary case)
    if (!is.null(primary_state_df)){
      # Get rows of primary_state_df matching primary cases
      infector_df <- subset(primary_state_df,case_id %in% primary_case_ids)
      # Calculate how many days between now and infector being isolated
      infector_df$d_until_iso <- infector_df$incubation_length +
        infector_df$isolation_delay -
        infector_df$days_infected
      # For each new case, determine the actual isolation delay
      iso_delay2 <- vapply(X=seq_along(iso_delay), FUN=function(ii){
        if (state_df[ii,]$is_traced){ # If traced, then check infector isolation time
          d_until_onset <- state_df$incubation_length[ii]
          d_until_pri_iso <- subset(infector_df,case_id==primary_case_ids[ii])$d_until_iso
          # If secondary case onset is *after* primary's isolation time, delay set to 0
          if (d_until_onset > d_until_pri_iso){
            return(0) # NB: return is to vapply
          } else{ # otherwise, delay is until primary case is isolated
            return(d_until_pri_iso-d_until_onset) # NB: return is to vapply
          }
        } else{ # If not traced, then use drawn iso_delay value
          return(iso_delay[ii]) # NB: return is to vapply
        }
      }, FUN.VALUE=999)
      # Update iso_delay vector to account for tracing
      iso_delay <- iso_delay2
    }
  }
  return(iso_delay)
}

#' Determine whether a case is being manually traced by public health
#'
#' Draws a uniform random number and compares to trace probability to determine
#' whether a case will be manually traced.
#'
#' The probability of a case being traced, \code{p_trace}, can either be a value that varies
#' with cluster size (if \code{sim_params$vary_trace} is \code{TRUE}) or a constant value equal
#' to \code{sim_params$p_trace} if \code{sim_params$vary_trace} is \code{FALSE}.
#'
#' When variable tracing based on cluster size is used, \code{n} is assumed to be the cluster size
#' because it is assumed that this function would be called be \code{create_state_df} for all
#' secondary infections together. The three possible variable tracing values are for cluster sizes of:
#' less than 10 cases, 11-24 cases and 25 or more cases. The vector \code{sim_params$p_trace_vary}
#' sets the trace probability for each of these three cluster sizes.
#'
#' @param n           Number of cases to consider
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, where
#'                    all of the information needed to determine the probability of a case being
#'                    traced (\code{p_trace}) are found within:
#'                    \itemize{
#'                    \item \code{sim_params$vary_trace}: Boolean to determine whether
#'                    \code{p_trace} varies based on cluster size or whether a single constant
#'                    value for \code{p_trace} is used.
#'                    \item \code{sim_params$p_trace_vary}: Defines the \code{p_trace} value for
#'                    each cluster size group. Ignored if \code{sim_params$vary_trace} is \code{FALSE}.
#'                    \item \code{sim_params$p_trace}: Constant \code{p_trace} value used at all times.
#'                    Ignored if \code{sim_params$vary_trace} is \code{TRUE}.
#'                    }
#' @return A boolean vector of length n for each case's manual tracing status
draw_traced_status <- function(n, sim_params){
  if (sim_params$vary_trace){
    if(n<10){
      p_trace<-sim_params$p_trace_vary[1]
    }
    else if(n<25){
      p_trace<-sim_params$p_trace_vary[2]
    } else{
      p_trace<-sim_params$p_trace_vary[3]
    }
    return(runif(n) < p_trace)
  } else{ # otherwise, use fixed tracing value
    return(runif(n) < sim_params$p_trace)
  }
}

#' Determine whether a case is using a contact tracing app
#'
#' Draws a uniform random number and compares to app usage probability to determine
#' whether a case is an app user. For a case to be traced via the app, the
#' secondary case and the primary/index case must both be app users.
#'
#' @param n           Number of cases to consider
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, here,
#'                    the \code{sim_params$p_trace_app} value is used as the probability of
#'                    contact tracing app usage.
#' @return A boolean vector of length n for each case's app usage status
draw_trace_app_user_status <- function(n, sim_params){
  return(runif(n) < sim_params$p_trace_app)
}

#' Helper function to check primary and secondary trace_app_user status
#'
#' For a case to be traced via the app, the secondary case and the
#' primary/index case must both be app users.
#'
#' @param state_df          \code{state_df} object for the new cases
#' @param primary_state_df  The \code{state_df} object for the index/primary cases.
#' @param primary_case_ids  A list of case_ids for the primary cases.
#' @return A boolean vector of length n for each case's app tracing status
draw_traced_by_app <- function(state_df, primary_state_df, primary_case_ids){
  infector_df <- subset(primary_state_df,case_id %in% primary_case_ids)
  infector_traced_by_app <- vapply(seq(1:nrow(state_df)),
                                   function(ii){
                                     subset(infector_df, case_id==primary_case_ids[ii])$is_trace_app_user
                                   },
                                   TRUE)
  return(state_df$is_trace_app_user & infector_traced_by_app)
}

#' Determine whether a case will comply with tracing app instructions
#'
#' Draws a uniform random number and compares to app compliance probability. This status
#' is determined (but never used) even for cases that aren't traced by the app.
#'
#' @param state_df    \code{state_df} object for the new cases
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, here,
#'                    the \code{sim_params$p_trace_app_comp} value is used as the probability of
#'                    contact tracing app compliance.
#' @return A boolean vector of length n for each case's app compliance status
draw_trace_app_compliance <- function(state_df, sim_params){
  n <- nrow(state_df)
  rng <- runif(n) < sim_params$p_trace_app_comp
  return(state_df$is_trace_app_user & rng)
}

#' Determine infection length of new cases
#'
#' Currently a constant number for all cases
#'
#' @param n           Number of cases required
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, here,
#'                    the \code{sim_params$infect_dur} value is used the infection length in days.
#' @return A vector of length n for infection duration (double)
draw_infection_length <- function(n, sim_params){
  return(rep(sim_params$infect_dur,n))
}

#' Assign physical distancing behaviour to new cases
#'
#' Distancing behaviour is defined as a number between 0 and 1, representing the
#' relative number of contacts the subject encounters (i.e. 0.6 means the subject
#' reduces contact to 60% of the non-distancing population.)
#'
#' Current implementation assume that either the entire population has the same
#' distancing behaviour or that the population can be divided into two groups:
#' one that is practicing distancing and one that is not. User provides parameters such
#' as the fraction of the general population that is distancing and the relative contact
#' rate of the distancing population. This relative contact rate can change once, at the
#' time point provided by the user.
#'
#' Note: The ratio of distancers to non-distancers of new cases is not the same as the
#' ratio amongst the general population (user-provided values). This is because new cases
#' do not arise from the general population but from the population that is in contact with
#' others. Therefore, the fraction of distancers being infected would generally be lower than
#' the fraction of distancers in the general population.
#'
#' @param n_cases     Number of cases required
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters.
#'                    Here, the sim_params$phys_dist_params object contains the information
#'                    needed. This list should have the following entries:
#'                    \itemize{
#'                    \item \code{pd_pop_frac}: A number for the fraction of the general population
#'                    that is practicing distancing. Set this to 0 or 1 to have everyone act the same.
#'                    \item \code{pd_contact_rate1}: The contact rate of the distancing population before
#'                    time \code{pd_change_t}. Set this to 1 to have the distancing population not begin
#'                    distancing until a later time.
#'                    \item \code{pd_contact_rate2}: The contact rate of the distancing population after
#'                    time \code{pd_change_t}. Set this to 1 to have the distancing population stop
#'                    distancing at a certain time.
#'                    \item \code{pd_change_t}: The simulation time where the distancing population's
#'                    contact rate changes. Set this to 0 or a very high number to have only one rate.
#'                    }
#' @param sim_status  \code{sim_status} object (a list) containing simulation state vector
#' @return A vector of length n for the contact_rate (double)
draw_contact_rate <- function(n_cases, sim_params, sim_status){
  pd_params <- sim_params$phys_dist_params
  # Allows for a different contact rate for the PD group before and after a change time
  if (sim_status$t < pd_params$pd_change_t){
    pd_contact_rate <- pd_params$pd_contact_rate1  # rate before change time
  } else {
    pd_contact_rate <- pd_params$pd_contact_rate2  # rate after change time
  }

  # One population only:
  if (pd_params$pd_pop_frac==1){ # Everyone is (or is not) distancing (i.e. just one population)
    contact_rate<-rep(pd_contact_rate, n_cases)
  }
  else{
    # For two groups, need to factor in how the groups encounter each other in order
    # to properly assign probabilities of **infecting** someone in a given group
    # Ratio of distancer to non-distancer: epsilon*lambda/(1-epsilon)
    # epsilion = fraction of distancing pop, lambda = contact rate of distancing pop.
    eps <- pd_params$pd_pop_frac
    lambda <- pd_contact_rate
    prob_pd <- eps*lambda / (eps*lambda + 1-eps)
    contact_rate <- sample(x = c(lambda, 1.0),
                        size = n_cases,
                        prob = c(prob_pd, 1-prob_pd),
                        replace=TRUE)
  }
  return(contact_rate)
}

#' Stage secondary infections when new cases are generated
#'
#' When a new case is generated, all the secondary infections are also generated and
#' loaded into the \code{state_df} row for the primary case. The infections are only
#' actually loaded into the simulation at the correct time, through the function
#' \code{generate_secondary_infections()} in generate_new_infections.R. To generate
#' secondary infections, first the number of potential secondary infections is drawn.
#' Then, generation intervals are generated for each potential secondary infection. Finally,
#' potential infections are rejected by user-defined criteria, such as being after the
#' isolation time for the primary case. Currently, only one type of secondary infection
#' algorithm is implemented (the one from Hellewell et al.). The generation intervals are
#' generated by a call to \code{draw_generation_interval}.
#'
#' @param state_df    \code{state_df} object for the newly generated cases
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters.
#'                    Specifically, the \code{sim_params$sec_infect_params} object contains
#'                    the parameters for generating secondary infections.
#'                    This list should have the following entries:
#'                    \itemize{
#'                    \item \code{type}: Defines the algorithm to be used. Currently, only
#'                    "Hellewell" is implemented. This uses a negative binomial distribution
#'                    to draw the number of potential secondary infections with a mean equal
#'                    to R0*contact_rate and a dispersion equal to the value set below. Generation
#'                    intervals are drawn but only potential infections with generation intervals
#'                    earlier than the primary case's time of isolation are kept.
#'                    \item \code{disp}: The dispersion value for the negative binomial distribution
#'                    used to determine the number of potential secondary infections.
#'                    }
#' @param sim_status  \code{sim_status} object (a list) containing simulation state vector
#' @param import      A boolean indicating whether these new cases are imported cases instead
#'                    of secondary infections. These cases may be required to self-isolate.
#'                    Defaults to \code{FALSE}.
#' @return A list containing
#' \itemize{
#'   \item \code{n} - A vector with number of accepted secondary infections for each case (int)
#'   \item \code{generation} - A list where each entry corresponds to each case and contains a
#'   vector of generation intervals of accepted secondary infections.
#'   \item \code{non_infects} - A list where each entry corresponds to each case and contains a
#'   vector of generation intervals of rejected secondary infections, for debugging.
#' }
#'
draw_sec_infects_df <- function(state_df, sim_params, sim_status, import=FALSE){
  n_cases = nrow(state_df)
  if (sim_params$sec_infect_params$type=='Hellewell'){
    # Following Hellewell et al, for each case, determine
    # number of sec. infections and generation interval of each sec. infection
    col_names <- c('n_infect', 'generation_int')
    n_cols <- length(col_names)
    # Determine which physical distancing group cases belong to
    contact_rate <- state_df$contact_rate
    # Determine number of secondary infections drawn from neg. binomial
    mean_infect <- sim_params$R0 * contact_rate
    disp_infect <- sim_params$sec_infect_params$disp
    n_infect <- rnbinom(n_cases, mu=mean_infect, size=disp_infect)
    # Determine generation interval of each secondary infection for each infector source case
    if (sim_params$generation_int_params$dist=='gamma'){ # BC case
      generation_int <- mapply(draw_generation_interval, # FUN to be called on each X
                           n_infect, # X to loop over
                           MoreArgs=list(sim_params=sim_params), # additional required input for FUN
                           SIMPLIFY = FALSE) # force return as list
    } else if (sim_params$generation_int_params$dist=='skew_norm'){
      generation_int <- mapply(draw_generation_interval, # FUN to be called on each X1, X2
                           n_infect,             # X1 to loop over
                           state_df$incubation_length, #X2 to loop over
                           MoreArgs = list(sim_params = sim_params),# additional required input for FUN
                           SIMPLIFY = FALSE) # force return as list
    }
    # Sort generation intervals in ascending order
    generation_int<-lapply(generation_int,sort)
    ## Determine the days that each case is contagious
    # First day based on whether it's imported (imported cases may be required to isolate)
    if (import){
      import_params <- sim_params$import_params
      first_day_contagious <- sample(x=import_params$iso_lengths,
                                     size=n_cases,
                                     prob=import_params$iso_p_group,
                                     replace=TRUE)
    } else {
      first_day_contagious <-rep(0,n_cases)
    }
    # Switches for each contagious scenario
    is_T_and_S <- state_df$is_traced * state_df$is_symptomatic # traced and symptomatic
    is_T_and_nS <- state_df$is_traced * (1-state_df$is_symptomatic) # traced and not sympt
    is_nT_and_S <- (1-state_df$is_traced) * state_df$is_symptomatic # not traced and sympt
    is_nT_and_nS <- (1-state_df$is_traced) * (1-state_df$is_symptomatic) # not traced and not sympt
    # Length of contagious period for each scenario
    T_and_S_time <- state_df$incubation_length + state_df$isolation_delay # delay for primary not yet isolated folded into isolation delay
    T_and_nS_time <- state_df$infection_length # no isolation if no symptoms
    nT_and_S_time <- state_df$incubation_length + state_df$isolation_delay # isolated some time after symptoms
    nT_and_nS_time <- state_df$infection_length # no isolation if no symptoms
    last_day_contagious <- is_T_and_S   * T_and_S_time   +
      is_T_and_nS  * T_and_nS_time  +
      is_nT_and_S  * nT_and_S_time  +
      is_nT_and_nS * nT_and_nS_time
    # Ensure no "last day contagious" is larger than infection_length days (would be removed as inactive)
    longer_than_inf_len <- last_day_contagious > state_df$infection_length
    last_day_contagious[longer_than_inf_len] <- state_df$infection_length[longer_than_inf_len]
    # Also, to avoid weird negative values in case somehow the isolated infector still caused
    # this case to happen, make all negative values be zero (would not allow any sec infects)
    neg_last_day <- last_day_contagious < 0
    last_day_contagious[neg_last_day] <- 0
    # Split the generation interval list.
    # Keep the valid infections for model, store the rest for record keeping
    generation_keep <- lapply(
      seq_along(generation_int),
      function(ii, generation, first_day, last_day){
        index_to_keep <- (generation[[ii]] > first_day[ii]) &
          (generation[[ii]] < last_day[ii])
        return(generation[[ii]][index_to_keep])
      },
      generation=generation_int,
      first_day=first_day_contagious,
      last_day=last_day_contagious
    )
    generation_reject <- lapply(
      seq_along(generation_int),
      function(ii, generation, first_day, last_day){
        index_to_keep <- (generation[[ii]] > first_day[ii]) &
          (generation[[ii]] < last_day[ii])
        return(generation[[ii]][!index_to_keep])
      },
      generation=generation_int,
      first_day=first_day_contagious,
      last_day=last_day_contagious
    )
    # Get number of actual infections
    n_infect <- sapply(generation_keep,length)
    # Will probably move this outside of the if-block when other methods added
    return(list(n=n_infect,generation=generation_keep,non_infects=generation_reject))
  }
}

#' Generate generation intervals for potential secondary infections
#'
#' Generation intervals for each case's potential secondary infections are drawn from
#' distribution defined in \code{sim_params$generation_int_params}
#'
#' @param n                   Number of potential secondary infections for this case
#' @param incubation_length   Incubation period of this case (used in skew normal distribution only)
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters, where
#'                    all of the information needed to describe the distribution is found within
#'                    \code{sim_params$generation_int_params}. Current distributions possible are:
#'                    \itemize{
#'                    \item \code{skew_norm}: Drawn from skew normal distribution with xi set
#'                    to the case's incubation period, and attributes "omega" and "alpha" given
#'                    in \code{sim_params$generation_int_params}. This is what Hellewell uses.
#'                    \item \code{gamma}: Drawn from a gamma distribution with attributes
#'                    "shape" and "rate" given in \code{sim_params$generation_int_params}.
#'                    }
#' @return A vector of length n for generation intervals of the case's potential secondary infections
draw_generation_interval <- function(n, incubation_length, sim_params){
  if (sim_params$generation_int_params$dist=='skew_norm'){
    sn_xi = incubation_length # case incubation period
    sn_omega = sim_params$generation_int_params$omega
    sn_alpha = sim_params$generation_int_params$alpha
    generation_ints <- sn::rsn(n, xi=sn_xi, omega=sn_omega, alpha=sn_alpha)
    return(as.numeric(generation_ints))
  } else if (sim_params$generation_int_params$dist=='gamma'){
    generation_ints <- rgamma(n, shape=sim_params$generation_int_params$shape, rate=sim_params$generation_int_params$rate)
    return(generation_ints)
  }
}
