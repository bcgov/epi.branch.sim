# generate_new_infections.R
# Contains functions that tells the model when to add new infections,
# either from secondary infection cases or imported cases.
# New or alternate infection methods should be added here.

#' Adds secondary infections to simulation
#'
#' Potential secondary infections are generated at the time the primary
#' case was generated. This function "loads" the infections into the simulation
#' at the correct time. It also modifies the \code{state_df} to account for the
#' loaded infection.
#'
#' @param state_df    \code{state_df} object for the simulation
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters
#' @return A list containing
#' \itemize{
#'   \item \code{state_df} - The updated \code{state_df} object for the simulation
#'   \item \code{sec_infection_sources} - A vector of \code{case_id}s causing the infections.
#'   If a case causes multiple secondary infections, the ID is repeated. For example,
#'   c(1,5,5,7) means case 1 caused 1 infection, case 5 caused 2 infections and case 7
#'   caused 1 infection. If no new secondary infections are generated in this time step,
#'   this value is \code{NULL}.
#' }
generate_secondary_infections <- function(state_df, sim.params){
  sec_infection_sources <- c()
  if (nrow(state_df)>0){
    # TODO: Can we replace this for-loop?
    for (row in 1:nrow(state_df)){
      serial_int <- state_df$serial_intervals[[row]]
      # Find which entries in serial_int are causing new secondary infections from this case
      sec_inf_ind <- state_df$days_infected[row] > serial_int
      n_sec_inf <- sum(sec_inf_ind)
      if (n_sec_inf < 1){
        next # no new infections caused by this case this step
      }
      # Add case_id to output list once for each new secondary infection case
      case_id <- state_df$case_id[row]
      sec_infection_sources <- c(sec_infection_sources,rep(case_id,n_sec_inf))
      # Decrease infection counter and remove from serial interval list
      state_df$n_sec_infects[row] <- state_df$n_sec_infects[row] - n_sec_inf
      state_df$serial_intervals[[row]] <- serial_int[!sec_inf_ind]
    }
    return(list(updated_state_df=state_df, sec_infection_sources=sec_infection_sources))
  }
  else{
    return(list(updated_state_df=state_df, sec_infection_sources=NULL))
  }
}

#' Adds imported infections (new cases) to simulation
#'
#' The daily risk of imported infections is defined from an external data file which
#' contains the mean number of infected cases arriving per day. This risk level is
#' loaded into the simulation as \code{sim.params$import_params}. This function "loads"
#' the infections into the simulation at integer time steps, drawing a number of new
#' infections as Poisson with mean equal to the daily risk.
#'
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters
#' @param sim.status  \code{sim.status} object (a list) containing simulation state vector
#' @return A list containing
#' \itemize{
#'   \item \code{n_imported_infections} - Integer with number of imported cases at this time step.
#'   \item \code{import_source} - String describing the source of these new cases. This information
#'   will be recorded in \code{state_df} and \code{record_df}. For now, this will always be "imported".
#' }
generate_imported_infections <- function(sim.params, sim.status){
  if (sim.params$import_params=="None"){
    n_import<-0
  }
  else{
    # Only import cases at integer timesteps
    if ((sim.status$t %% 1) == 0){
      risk_index <- sim.status$t
      risk<-sim.params$import_params[risk_index]
      n_import <- rpois(1,risk)
    }
    else {
      n_import<-0
    }
  }
  return(list(n_imported_infections=n_import, import_source='imported'))
}
