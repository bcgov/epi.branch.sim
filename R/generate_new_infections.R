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
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters
#' @return A list containing
#' \itemize{
#'   \item \code{state_df} - The updated \code{state_df} object for the simulation
#'   \item \code{sec_infection_sources} - A vector of \code{case_id}s causing the infections.
#'   If a case causes multiple secondary infections, the ID is repeated. For example,
#'   c(1,5,5,7) means case 1 caused 1 infection, case 5 caused 2 infections and case 7
#'   caused 1 infection. If no new secondary infections are generated in this time step,
#'   this value is \code{NULL}.
#' }
generate_secondary_infections <- function(state_df, sim_params){
  sec_infection_sources <- c()
  if (nrow(state_df)>0){
    # TODO: Can we replace this for-loop?
    for (row in 1:nrow(state_df)){
      generation_int <- state_df$generation_intervals[[row]]
      # Find which entries in generation_int are causing new secondary infections from this case
      sec_inf_ind <- state_df$days_infected[row] > generation_int
      n_sec_inf <- sum(sec_inf_ind)
      if (n_sec_inf < 1){
        next # no new infections caused by this case this step
      }
      # Add case_id to output list once for each new secondary infection case
      case_id <- state_df$case_id[row]
      sec_infection_sources <- c(sec_infection_sources,rep(case_id,n_sec_inf))
      # Decrease infection counter and remove from generation interval list
      state_df$n_sec_infects[row] <- state_df$n_sec_infects[row] - n_sec_inf
      state_df$generation_intervals[[row]] <- generation_int[!sec_inf_ind]
    }
    return(list(updated_state_df=state_df, sec_infection_sources=sec_infection_sources))
  }
  else{
    return(list(updated_state_df=state_df, sec_infection_sources=NULL))
  }
}

#' Adds imported infections to simulation
#'
#' The number of imported infections may be a deterministic or stochastic value.
#'
#' The list object \code{sim_params$import_params} determines how imported cases are
#' added to the simulation. The following import types are supported, which is chosen by
#' \code{sim_params$import_params$type}:
#' \describe{
#'  \item{\code{type}="None" or "none"}{No cases are imported.}
#'
#'  \item{\code{type}="constant"}{Each day, the number of imported cases is equal to a constant value, set
#'  by \code{sim_params$import_params$rate}.}
#'
#'  \item{\code{type}="constant_two_phase"}{Same as the above, but with two possible constant values. For
#'  times prior to \code{sim_params$import_params$delay}, the rate is \code{sim_params$import_params$rate1}.
#'  After this time, the rate is \code{sim_params$import_params$rate2}.}
#'
#'  \item{\code{type}="daily_constant"}{Each day, the number of imported cases is given by a vector defined by
#'  \code{sim_params$import_params$rate}, which should have length equal to the number of days in the
#'  simulation.}
#'
#'  \item{\code{type}="poisson"}{Each day, the number of imported cases is drawn from a Poisson distribution
#'  with mean set by \code{sim_params$import_params$rate}.}
#'
#'  \item{\code{type}="daily_risk"}{Each day, the number of imported cases is drawn from a Poisson
#'  distribution with mean equal to the daily risk. The daily risk is given by a vector defined by
#'  \code{sim_params$import_params$risk}, which should have length equal to the number of days in the
#'  simulation.}
#' }
#'
#' @param sim_params  \code{sim_params} object (a list) containing simulation parameters
#' @param sim_status  \code{sim_status} object (a list) containing simulation state vector
#' @return A list containing
#' \itemize{
#'   \item \code{n_imported_infections} - Integer with number of imported cases at this time step.
#'   \item \code{import_source} - String describing the source of these new cases. This information
#'   will be recorded in \code{state_df} and \code{record_df}. For now, this will always be "imported".
#' }
generate_imported_infections <- function(sim_params, sim_status){
  import_params <- sim_params$import_params
  if (import_params$type=="None" | import_params$type=='none'){
    n_import<-0
  }
  else if (import_params$type=='constant'){
    n_import <- import_params$rate
  }
  else if (import_params$type=='constant_two_phase'){
    if (sim_status$t < import_params$delay){
      n_import <- import_params$rate1
    } else{
      n_import <- import_params$rate2
    }
  }
  else if (import_params$type=='daily_constant'){
    rate_index <- as.integer(round(sim_status$t))
    n_import <- import_params$rate[rate_index]
  }
  else if (import_params$type=='poisson'){
    n_import <- rpois(1,import_params$rate)
  }
  else if (import_params$type=='daily_risk'){
    risk_index <- as.integer(round(sim_status$t))
    risk<-import_params$risk[risk_index]
    n_import <- rpois(1,risk)
  }
  return(list(n_imported_infections=n_import, import_source='imported'))
}
