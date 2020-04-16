# These functions are used to initialize simulation control objects

#' Creates \code{sim.params} object for global simulation parameters
#'
#' This object is a list containing global simulation parameters such as the disease
#' parameters, either directly defined or through a distribution with user-provided
#' values. This object also contains some parameters that control the simulation itself.
#' Many of these values are used in other functions, the major ones are linked here.
#'
#' @param R0            Reproductive number of disease. Used in \code{\link{draw_sec_infects_df}} to
#'                      determine the number of potential secondary infections
#'                      from a primary case.
#' @param infect.dur    Duration of infection (in days). This constant value determines when
#'                      a case is no longer active. It's also the length of the infectious
#'                      period, unless the case is isolated. See \code{\link{draw_infection_length}}
#' @param vary.trace    Boolean to determine whether contact tracing is a constant or variable
#'                      based on the number of new cases in a cluster.
#'                      See \code{\link{draw_traced_status}}
#' @param p.trace,cluster_breaks,cluster_p.trace Parameters determining contact tracing
#'                                               probability. See \code{\link{draw_traced_status}}
#' @param p.symp        Probability that a case is symptomatic.
#'                      Used in \code{\link{draw_symptomatic_status}}.
#' @param dt            Size of the simulation time step (in days).
#' @param incub_params  Parameters for distribution of incubation length.
#'                      See \code{\link{draw_incubation_period}}.
#' @param serial_int_params  Parameters for distribution of serial intervals.
#'                           See \code{\link{draw_serial_interval}}.
#' @param iso_delay_params   Parameters for distribution of delay from symptom onset to isolation.
#'                           See \code{\link{draw_isolation_delay_period}}.
#' @param sec_infect_params  Parameters to model secondary infections.
#'                           See \code{\link{draw_sec_infects_df}}
#' @param import_params      Parameters to model imported infections.
#'                           See \code{\link{generate_imported_infections}}
#' @param social_dist_params Parameters to model social/physical distancing or
#'                           other interventions that reduce the potential to
#'                           cause secondary infections.
#'                           See \code{\link{draw_sd_factor}}
#' @return A list containing all of the above entries
initialize_sim_params <- function(R0, infect.dur, vary.trace, p.trace, cluster_breaks,
                                  cluster_p.trace, p.symp, dt,
                                  incub_params, serial_int_params,
                                  iso_delay_params, sec_infect_params,
                                  import_params, social_dist_params){
  sim.params <- list(
    R0=R0,                      # used to parameterize probability of infection
    infect.dur=infect.dur,      # typical infection duration
    vary.trace=vary.trace,      # Boolean on whether using variable tracing model
    p.trace=p.trace,            # probability of being traced (fixed, ignored if vary.trace is TRUE)
    cluster_breaks=cluster_breaks,   # breakpoints for cluster sizes (ignored if vary.trace is FALSE)
    cluster_p.trace=cluster_p.trace, # p.trace for each cluster size group (ignored if vary.trace is FALSE)
    p.symp=p.symp,              # probability of symptomatic case
    dt=dt,                      # size of time step
    incub_params=incub_params,           # parameters for distribution of incubation length
    serial_int_params=serial_int_params, # parameters for distribution of serial interval
    iso_delay_params=iso_delay_params,   # parameters for distribution of isolation delay length
    sec_infect_params=sec_infect_params,  # parameters to model secondary infections
    # Daily import hazards
    import_params=import_params,
    # Social distancing parameters
    social_dist_params=social_dist_params
  )
  return(sim.params)
}

#' Creates \code{sim.status} object, the simulation state vector
#'
#' This object is a list containing the simulation state vector. This contains values that
#' change with every time step. Note that \code{sim.params} is another important object
#' containing quantities that remain constant throughout the simulation. Currently, the
#' only two items in \code{sim.status} are the time and the case_id of the last case
#' infected (used to generate the next case_id). This function creates the \code{sim.status}
#' object at the beginning of the simulation and only needs to be called once.
#'
#' @param start_time      The start time (in days) for the simulation
#' @param initial_n_cases The number of cases at the start of the simulation.
#' @return A list containing the simulation state vector quantities, which includes:
#' \itemize{
#'  \item \code{t} Current time (days)
#'  \item \code{last_case_id} The last case_id generated (integer)
#' }
initialize_sim_status <- function(start_time, initial_n_cases){
  sim.status <- list(
    t=start_time,
    last_case_id=initial_n_cases
  )
  return(sim.status)
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
#' \code{is_symptomatic}, \code{days_infected}, \code{incubation_length},
#' \code{isolation_delay}, \code{infection_length}, \code{sd_factor},
#' \code{n_sec_infects}, and \code{serial_intervals}. The meaning of these columns
#' will be documented elsewhere. Most of the columns are constant and are initialized
#' in this function via calls to related functions in draw_distributions.R. All time lengths
#' are in days and represent a relative time, not absolute times. The
#' \code{status}, \code{n_sec_infects} and \code{serial_intervals} columns are updated
#' over the course of the simulation as cases progress through infection states. The
#' \code{n_sec_infects} counter and corresponding \code{serial_intervals} list have its
#' elements removed as the secondary infections are loaded into the simulation. As cases
#' become isolated or inactive and no longer cause infections, they are removed from this
#' dataframe.
#'
#' A companion dataframe, \code{rec_df}, is meant to be a complete record of all cases
#' that occurred in the simulation. This object is created by \code{\link{create_record_df}}
#' and tracks all times as the absolute amount of time passed since the simulation start.
#'
#' @param n_cases      The number of cases to initialize in this dataframe
#' @param sim.params   The \code{sim.params} object containing disease and simulation parameters
#' @param sim.status   The \code{sim.status} state vector
#' @param initialize   A boolean indicating whether these cases are the first cases of the
#'                     simulation. Defaults to \code{FALSE}.
#' @param import       A boolean indicating whether these new cases are imported cases instead
#'                     of secondary infections. Defaults to \code{FALSE}.
#' @return A \code{state_df} dataframe object
create_state_df <- function(n_cases, sim.params, sim.status, initialize=FALSE, import=FALSE){
  # List of columns in state df
  col_names <- c("case_id", "status", "is_traced", "is_symptomatic", "days_infected",
                 "incubation_length", "isolation_delay", "infection_length", "sd_factor",
                 "n_sec_infects", "serial_intervals") #
  n_cols <- length(col_names)
  # Create data frame
  state_df <- data.frame(matrix(nrow=n_cases,ncol=n_cols, dimnames=list(NULL,col_names)))
  # Fill in start values
  state_df$case_id <- sim.status$last_case_id + 1:n_cases
  state_df$status <- rep("incubation", n_cases)
  if (initialize){
    state_df$is_traced <- rep(FALSE, n_cases)
  } else{
    state_df$is_traced <- draw_traced_status(n_cases,sim.params)
  }
  state_df$is_symptomatic <- draw_symptomatic_status(n_cases,sim.params)
  state_df$days_infected <- rep(0, n_cases)
  state_df$incubation_length <- draw_incubation_period(n_cases,sim.params)
  state_df$isolation_delay <- draw_isolation_delay_period(state_df,sim.params)
  state_df$infection_length <- draw_infection_length(n_cases, sim.params)
  sec_infect_out <- draw_sec_infects_df(state_df,sim.params, sim.status, initialize=initialize, import=import)
  state_df$sd_factor <- sec_infect_out$sd_factor
  state_df$n_sec_infects<- sec_infect_out$n
  state_df$serial_intervals<- sec_infect_out$serial
  state_df$non_infect_serials <- sec_infect_out$non_infects
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
#' \code{is_symptomatic}, \code{d.incub}, \code{d.iso_delay}, \code{d.infection},
#' \code{sd_factor}, \code{n_sec_infects}, \code{d.serial_ints}, \code{t.inf},
#' \code{t.symp}, \code{t.iso}, \code{t.inact}, \code{cases_inf}, \code{s.status},
#' and \code{non_infect_serials}. The meaning of these columns will be documented elsewhere.
#' Columns \code{case_id} through \code{t.inf} as well as column
#' \code{non_infect_serials} are constant and are initialized when this function
#' is called, using the provided \code{state_df} object. The columns \code{t.symp},
#' \code{t.iso} and \code{t.inact} are initialized as \code{NA} and replaced with the actual
#' simulation time where these events happen. If they do not happen (e.g. an asymptomatic
#' case will not have a \code{t.symp} or \code{t.iso}), then the value remains as \code{NA}.
#' The column \code{cases_inf} is meant for a list of case_id of all secondary infections
#' caused by the given case but is not yet implemented. The column \code{s.status} is updated
#' by \code{\link{step_simulation}} at the appropriate time in the simulation.
#'
#' Note that all time lengths are given in days and values in columns beginning with
#' \code{d.} are a period of days since the relevant reference time while columns
#' beginning with \code{t.} are absolute time, measured in days since the simulation began.
#'
#' @param state_df     \code{state_df} object for the cases to be generated
#' @param sim.status   The \code{sim.status} state vector
#' @param initialize   A boolean indicating whether these cases are the first cases of the
#'                     simulation. Defaults to \code{FALSE}.
#' @param infection_source  A string or vector of case_ids representing the source of each
#'                          case's infection. Defaults to \code{NULL}, however, this
#'                          parameter cannot be NULL if \code{initialize} is \code{FALSE}.
#' @return A \code{rec_df} dataframe object
create_record_df <- function(state_df, sim.status, initialize=FALSE, infection_source=NULL){
  # List of columns to record
  col_names <- c("case_id", "source", "is.traced", "is.symptomatic", "d.incub", "d.iso_delay",
                 "d.infection", "sd_factor", "n.sec_infects", "d.serial_ints",
                 "t.inf", "t.symp", "t.iso", "t.inact",
                 "cases_inf", "s.status", "non_infect_serials")
  n_cols <- length(col_names)
  # Get number of rows
  n_rows <- nrow(state_df)
  # Create data frame
  rec_df <- data.frame(matrix(nrow=n_rows,ncol=n_cols, dimnames=list(NULL,col_names)))
  # Populate data frame
  rec_df$case_id <- state_df$case_id
  if (initialize){
    rec_df$source <- rep("initial", n_rows) # initial cases
  } else{
    rec_df$source <- infection_source # vector of source case_ids or import source provided
  }
  rec_df$is.traced <- state_df$is_traced
  rec_df$is.symptomatic <- state_df$is.symptomatic
  rec_df$d.incub <- state_df$incubation_length
  rec_df$d.iso_delay <- state_df$isolation_delay
  rec_df$d.infection <- state_df$infection_length
  rec_df$sd_factor <- state_df$sd_factor
  rec_df$n.sec_infects <- state_df$n_sec_infects
  rec_df$d.serial_ints <- state_df$serial_intervals
  rec_df$t.inf <- rep(sim.status$t, n_rows)
  rec_df$s.status <- state_df$status
  rec_df$non_infect_serials <- state_df$non_infect_serials

  return(rec_df)
}

