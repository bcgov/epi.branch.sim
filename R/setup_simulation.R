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

#' Creates \code{sim_params} object for global simulation parameters
#'
#' This object is a list containing global simulation parameters such as the disease
#' parameters, either directly defined or through a distribution with user-provided
#' values. This object also contains some parameters that control the simulation itself.
#' Many of these values are used in other functions, the major ones are linked here.
#'
#' @param R0            Reproductive number of disease. Used in \code{\link{draw_sec_infects_df}} to
#'                      determine the number of potential secondary infections
#'                      from a primary case.
#' @param infect_dur    Duration of infection (in days). This constant value determines when
#'                      a case is no longer active. It's also the length of the infectious
#'                      period, unless the case is isolated. See \code{\link{draw_infection_length}}
#' @param vary_trace    Boolean to determine whether contact tracing is a constant or variable
#'                      based on the number of new cases in a cluster.
#'                      See \code{\link{draw_traced_status}}
#' @param p_trace       Probability that a case will be manually traced.
#' @param p_trace_app   Probability that a case will be a contact tracing app user.
#' @param p_trace_app_comp Probability that a case would comply with contact tracing app recommendations.
#' @param p_symp        Probability that a case is symptomatic.
#'                      Used in \code{\link{draw_symptomatic_status}}.
#' @param dt            Size of the simulation time step (in days).
#' @param incub_params  Parameters for distribution of incubation length.
#'                      See \code{\link{draw_incubation_period}}.
#' @param generation_int_params  Parameters for distribution of generation intervals.
#'                           See \code{\link{draw_generation_interval}}.
#' @param iso_delay_params   Parameters for distribution of delay to isolation.
#'                           See \code{\link{draw_isolation_delay_period}}.
#' @param sec_infect_params  Parameters to model secondary infections.
#'                           See \code{\link{draw_sec_infects_df}}
#' @param import_params      Parameters to model imported infections.
#'                           See \code{\link{generate_imported_infections}}
#' @param phys_dist_params   Parameters to model physical (social) distancing or
#'                           other interventions that reduce the potential to
#'                           cause secondary infections.
#'                           See \code{\link{draw_contact_rate}}
#' @return A list containing all of the above entries
#' @export
initialize_sim_params <- function(R0, infect_dur, vary_trace, p_trace,
                                  p_trace_app, p_trace_app_comp, p_symp, dt,
                                  incub_params, generation_int_params,
                                  iso_delay_params, sec_infect_params,
                                  import_params, phys_dist_params){
  sim_params <- list(
    R0=R0,                      # used to parameterize probability of infection
    infect_dur=infect_dur,      # typical infection duration
    vary_trace=vary_trace,      # Boolean on whether using variable tracing model
    p_trace=p_trace,            # probability of being traced (fixed, ignored if vary_trace is TRUE)
    p_trace_vary=c(1,0.8,0),    # probability of being traced for varying cluster sizes (ignored if vary_trace is FALSE)
    p_trace_app=p_trace_app,    # probability that a person is an app user
    p_trace_app_comp=p_trace_app_comp,    # probability that a person listens to app
    p_symp=p_symp,              # probability of symptomatic case
    dt=dt,                      # size of time step
    incub_params=incub_params,           # parameters for distribution of incubation length
    generation_int_params=generation_int_params, # parameters for distribution of generation interval
    iso_delay_params=iso_delay_params,   # parameters for distribution of isolation delay length
    sec_infect_params=sec_infect_params,  # parameters to model secondary infections
    # Daily import hazards
    import_params=import_params,
    # Physical distancing parameters
    phys_dist_params=phys_dist_params
  )
  return(sim_params)
}

#' Creates \code{sim_status} object, the simulation state vector
#'
#' This object is a list containing the simulation state vector. This contains values that
#' change with every time step. Note that \code{sim_params} is another important object
#' containing quantities that remain constant throughout the simulation. Currently, the
#' only two items in \code{sim_status} are the time and the case_id of the last case
#' infected (used to generate the next case_id). This function creates the \code{sim_status}
#' object at the beginning of the simulation and only needs to be called once.
#'
#' @param start_time      The start time (in days) for the simulation
#' @param start_n_cases   The number of cases at the start of the simulation.
#' @return A list containing the simulation state vector quantities, which includes:
#' \itemize{
#'  \item \code{t} Current time (days)
#'  \item \code{last_case_id} The last case_id generated (integer)
#' }
#' @export
initialize_sim_status <- function(start_time, start_n_cases){
  sim_status <- list(
    t=start_time,
    last_case_id=start_n_cases
  )
  return(sim_status)
}

#' Creates a \code{state_df} object
#'
#' This function creates a \code{state_df} dataframe object, which is used to track
#' active cases in the simulation.
#'
#' The simulation's main \code{state_df} dataframe is created with this
#' function at the beginning of the simulation. When new cases are added (secondary or
#' imported infections), a \code{state_df} dataframe is created for these new cases
#' by \code{\link{generate_secondary_infections}} and \code{\link{generate_imported_infections}}.
#' The new \code{state_df} object is merged with the main \code{state_df} within
#' \code{\link{step_simulation}}.
#'
#' The columns of this dataframe are: \code{case_id}, \code{status}, \code{is_traced},
#' \code{is_trace_app_user}, \code{is_trace_app_comply}, \code{is_traced_by_app},
#' \code{is_symptomatic}, \code{days_infected}, \code{incubation_length},
#' \code{isolation_delay}, \code{infection_length}, \code{contact_rate},
#' \code{n_sec_infects}, and \code{generation_intervals}. The meaning of these columns
#' will be documented elsewhere. Most of the columns are constant and are initialized
#' in this function via calls to related functions in draw_distributions.R. All time lengths
#' are in days and represent a relative time, not absolute times. The
#' \code{status}, \code{n_sec_infects} and \code{generation_intervals} columns are updated
#' over the course of the simulation as cases progress through infection states. The
#' \code{n_sec_infects} counter and corresponding \code{generation_intervals} list have its
#' elements removed as the secondary infections are loaded into the simulation. As cases
#' become isolated or inactive and no longer cause infections, they are removed from this
#' dataframe.
#'
#' A companion dataframe, \code{rec_df}, is meant to be a complete record of all cases
#' that occurred in the simulation. This object is created by \code{\link{create_record_df}}
#' and tracks all times as the absolute amount of time passed since the simulation start.
#'
#' @param n_cases      The number of cases to initialize in this dataframe
#' @param sim_params   The \code{sim_params} object containing disease and simulation parameters
#' @param sim_status   The \code{sim_status} state vector
#' @param initialize   A boolean indicating whether these cases are the first cases of the
#'                     simulation. Defaults to \code{FALSE}. Initial cases are never traced.
#' @param import       A boolean indicating whether these new cases are imported cases instead
#'                     of secondary infections. Defaults to \code{FALSE}. Imported cases are never
#'                     app users and are never traced.
#' @param primary_state_df  When generating a \code{state_df} object for secondary infections, this
#'                          parameter should be the \code{state_df} object for the index/primary cases.
#'                          This should be \code{NULL} when there is no primary case. Defaults to \code{NULL}.
#' @param primary_case_ids  When generating a \code{state_df} object for secondary infections, this
#'                          parameter should be a list of case_ids for the primary cases.
#'                          This should be \code{NULL} when there is no primary case. Defaults to \code{NULL}.
#' @return A \code{state_df} dataframe object
#' @export
create_state_df <- function(n_cases, sim_params, sim_status,
                            initialize=FALSE, import=FALSE,
                            primary_state_df=NULL, primary_case_ids=NULL){
  # List of columns in state df
  col_names <- c("case_id", "status", "is_traced", "is_trace_app_user", "is_trace_app_comply",
                 "is_traced_by_app", "is_symptomatic", "days_infected", "incubation_length",
                 "isolation_delay",  "infection_length", "contact_rate", "n_sec_infects",
                 "generation_intervals") #
  n_cols <- length(col_names)
  # Create data frame
  state_df <- data.frame(matrix(nrow=n_cases,ncol=n_cols, dimnames=list(NULL,col_names)))

  # If no cases being added, then return empty data frame
  # This should only happen when initializing
  if (n_cases==0){
    return(state_df)
  }

  # Fill in start values
  if (initialize){
    state_df$case_id <- 1:n_cases # special case for starting out
  } else{
    state_df$case_id <- sim_status$last_case_id + 1:n_cases
  }
  state_df$status <- rep("incubation", n_cases)
  if (initialize || import){
    state_df$is_traced <- rep(FALSE, n_cases)
  } else{
    state_df$is_traced <- draw_traced_status(n_cases,sim_params)
  }
  if (import){ # imported cases are not app users
    state_df$is_trace_app_user <- rep(FALSE, n_cases)
    state_df$is_trace_app_comply <- rep(FALSE, n_cases)
  } else{
    state_df$is_trace_app_user <- draw_trace_app_user_status(n_cases, sim_params)
    state_df$is_trace_app_comply <- draw_trace_app_compliance(state_df, sim_params)
  }
  if (!is.null(primary_state_df)){
    state_df$is_traced_by_app <- draw_traced_by_app(state_df, primary_state_df, primary_case_ids)
  } else {
    state_df$is_traced_by_app <- rep(FALSE, n_cases)
  }
  state_df$is_symptomatic <- draw_symptomatic_status(n_cases,sim_params)
  state_df$contact_rate <- draw_contact_rate(n_cases, sim_params, sim_status)
  state_df$days_infected <- rep(0, n_cases)
  state_df$incubation_length <- draw_incubation_period(n_cases,sim_params)
  state_df$isolation_delay <- draw_isolation_delay_period(state_df,sim_params,
                                                          primary_state_df=primary_state_df,
                                                          primary_case_ids=primary_case_ids)
  state_df$infection_length <- draw_infection_length(n_cases, sim_params)
  sec_infect_out <- draw_sec_infects_df(state_df,sim_params, sim_status, import=import)
  state_df$n_sec_infects<- sec_infect_out$n
  state_df$generation_intervals<- sec_infect_out$generation
  state_df$non_infect_generations <- sec_infect_out$non_infects
  # Return
  return(state_df)
}

#' Creates a \code{rec_df} object
#'
#' This function creates a \code{rec_df} dataframe object, which is used to record all
#' cases that ever get created within a simulation.
#'
#' The simulation's main \code{rec_df} dataframe is created with this
#' function at the beginning of the simulation. When new cases are added (secondary or
#' imported infections), a \code{rec_df} dataframe is created for these new cases
#' by \code{\link{generate_secondary_infections}} and \code{\link{generate_imported_infections}}.
#' The new \code{rec_df} object is merged with the main \code{rec_df} within
#' \code{\link{step_simulation}}.
#'
#' The columns of this dataframe are: \code{case_id}, \code{source}, \code{is_traced},
#' \code{is_trace_app_user}, \code{is_traced_by_app}, \code{is_trace_app_comply},
#' \code{is_symptomatic}, \code{d_incub}, \code{d_iso_delay}, \code{d_infection},
#' \code{contact_rate}, \code{n_sec_infects}, \code{d_generation_ints}, \code{t_inf},
#' \code{t_symp}, \code{t_iso}, \code{t_inact}, \code{cases_inf}, \code{s_status},
#' and \code{non_infect_generations}. The meaning of these columns will be documented elsewhere.
#' Columns \code{case_id} through \code{t_inf} as well as column
#' \code{non_infect_generations} are constant and are initialized when this function
#' is called, using the provided \code{state_df} object. The columns \code{t_symp},
#' \code{t_iso} and \code{t_inact} are initialized as \code{NA} and replaced with the actual
#' simulation time where these events happen. If they do not happen (e.g. an asymptomatic
#' case will not have a \code{t_symp} or \code{t_iso}), then the value remains as \code{NA}.
#' The column \code{cases_inf} is meant for a list of case_id of all secondary infections
#' caused by the given case but is not yet implemented. The column \code{s_status} is updated
#' by \code{\link{step_simulation}} at the appropriate time in the simulation.
#'
#' Note that all time lengths are given in days and values in columns beginning with
#' \code{d_} are a period of days since the relevant reference time while columns
#' beginning with \code{t_} are absolute time, measured in days since the simulation began.
#'
#' @param state_df     \code{state_df} object for the cases to be generated
#' @param sim_status   The \code{sim_status} state vector
#' @param initialize   A boolean indicating whether these cases are the first cases of the
#'                     simulation. Defaults to \code{FALSE}.
#' @param infection_source  A string or vector of case_ids representing the source of each
#'                          case's infection. Defaults to \code{NULL}, however, this
#'                          parameter cannot be NULL if \code{initialize} is \code{FALSE}.
#' @return A \code{rec_df} dataframe object
#' @export
create_record_df <- function(state_df, sim_status, initialize=FALSE, infection_source=NULL){
  # List of columns to record
  col_names <- c("case_id", "source", "is_traced", "is_trace_app_user", "is_traced_by_app",
                 "is_trace_app_comply", "is_symptomatic", "d_incub", "d_iso_delay",
                 "d_infection", "contact_rate", "n_sec_infects", "d_generation_ints",
                 "t_inf", "t_symp", "t_iso", "t_inact",
                 "cases_inf", "s_status", "non_infect_generations")
  n_cols <- length(col_names)
  # Get number of rows
  n_rows <- nrow(state_df)
  # Create data frame
  rec_df <- data.frame(matrix(nrow=n_rows,ncol=n_cols, dimnames=list(NULL,col_names)))

  # If no cases being added, then return empty data frame
  # This should only happen when initializing
  if (n_rows==0){
    return(rec_df)
  }

  # Populate data frame
  rec_df$case_id <- state_df$case_id
  if (initialize){
    rec_df$source <- rep("initial", n_rows) # initial cases
  } else{
    rec_df$source <- infection_source # vector of source case_ids or import source provided
  }
  rec_df$is_traced <- state_df$is_traced
  rec_df$is_trace_app_user <- state_df$is_trace_app_user
  rec_df$is_trace_app_comply <- state_df$is_trace_app_comply
  rec_df$is_traced_by_app <- state_df$is_traced_by_app
  rec_df$is_symptomatic <- state_df$is_symptomatic
  rec_df$d_incub <- state_df$incubation_length
  rec_df$d_iso_delay <- state_df$isolation_delay
  rec_df$d_infection <- state_df$infection_length
  rec_df$contact_rate <- state_df$contact_rate
  rec_df$n_sec_infects <- state_df$n_sec_infects
  rec_df$d_generation_ints <- state_df$generation_intervals
  #rec_df$sec_infects <- state_df$sec_infects
  rec_df$t_inf <- rep(sim_status$t, n_rows)
  rec_df$s_status <- state_df$status
  rec_df$non_infect_generations <- state_df$non_infect_generations

  return(rec_df)
}

