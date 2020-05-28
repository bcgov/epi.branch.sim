#' Takes the simulation one dt timestep forward
#'
#' Use the outputs of this function as input for the next step.
#'
#' @param sim.status  \code{sim.status} list containing simulation state vector
#' @param state_df    \code{state_df} object for the simulation
#' @param rec_df      \code{record_df} object for the simulation
#' @param sim.params  \code{sim.params} object (a list) containing simulation parameters
#' @return A list containing
#' \itemize{
#'   \item \code{status}   - The updated \code{sim.status} object
#'   \item \code{state}    - The updated \code{state_df} object
#'   \item \code{record}   - The updated \code{record_df} object
#'   \item \code{new_sec_cases} - Number of new secondary infections generated in this time step
#'   \item \code{new_imp_cases} - Number of new imported infectoins generated in this time step
#' }
step_simulation <- function(sim.status, state_df, rec_df, sim.params){
  # Increment time
  sim.status$t <- sim.status$t + sim.params$dt
  # Update timestamp of current infections
  state_df$days_infected <- state_df$days_infected + sim.params$dt

  # Determine which cases will cause infections this step
  gen_sec_infect_out <- generate_secondary_infections(state_df, sim.params)
  state_df <- gen_sec_infect_out$updated_state_df
  sec_infs_source <- gen_sec_infect_out$sec_infection_sources
  n_sec_infections <- length(sec_infs_source)
  # If new infections, create a temporary state and rec dataframes for the new cases
  # Also update rec_df for the origin cases and sim.status for case_id info
  if (n_sec_infections > 0){
    # Create new state_df for new cases
    sec_cases_state_df <- create_state_df(n_sec_infections,
                                          sim.params,
                                          sim.status,
                                          primary_state_df = state_df,
                                          primary_case_ids = sec_infs_source)
    # Update last_case_id
    sim.status$last_case_id <- max(sec_cases_state_df$case_id)
    # Create new rec_df for new cases
    sec_cases_rec_df <- create_record_df(sec_cases_state_df, sim.status, infection_source=sec_infs_source)
    # TODO: Update rec_df of the origin cases (so no tracking of this yet)
  }

  # Determine whether there will be any imported cases
  # Only import cases at whole number timesteps
  if ( (sim.status$t%%1) < 1e-3 || 1-(sim.status$t%%1) < 1e-3 ){
    gen_imp_infect_out <- generate_imported_infections(sim.params, sim.status)
    n_imp_infections <- gen_imp_infect_out$n_imported_infections
    imp_source <- gen_imp_infect_out$import_source
  } else{
    n_imp_infections<-0
  }
  # If there are imported cases, create temporary state and rec dataframes for the new cases
  # Also update sim.status for case_id info
  if (n_imp_infections > 0){
    # Create new state_df for new imported cases
    imp_cases_state_df <- create_state_df(n_imp_infections,
                                          sim.params,
                                          sim.status,
                                          import=TRUE)
    # Update last_case_id
    sim.status$last_case_id <- max(imp_cases_state_df$case_id)
    # Create new rec_df for new imported cases
    imp_cases_rec_df <- create_record_df(imp_cases_state_df, sim.status, infection_source=imp_source)
  }

  # First, identify and mark all cases where infection has ended, regardless of status
  cases_to_inact <- state_df$days_infected > state_df$infection_length
  n_cases_to_inact <- sum(cases_to_inact)
  if (n_cases_to_inact > 0){
    # Update state_df
    state_df$status[cases_to_inact] <- "inactive"
    # Update rec_df
    ids_to_inact <- state_df$case_id[cases_to_inact]
    rec_df$t.inact[rec_df$case_id %in% ids_to_inact] <- sim.status$t
    rec_df$s.status[rec_df$case_id %in% ids_to_inact] <- "inactive"
  }

  # Identify all remaining cases that are eligible for advancement to next stage
  # Doing it now prevents checking advanced cases for another advancement
  cases_adv_inc <- (state_df$status=="incubation") &
    (state_df$days_infected > state_df$incubation_length)
  cases_adv_symp <- (state_df$status=="symptomatic") &
    (state_df$days_infected > (state_df$incubation_length + state_df$isolation_delay))

  # Advance current infections from incubation to next stage
  # Goes to symptomatic or asymptomatic based on $is_symptomatic
  # Symptomatic cases that are traced goes right to isolation
  n_cases_adv_inc <- sum(cases_adv_inc)
  if (n_cases_adv_inc > 0){
    # Check whether advancing right to isolation (i.e. delay short enough)
    past_iso_delay <- state_df$days_infected > (state_df$incubation_length + state_df$isolation_delay)
    # Update state_df
    to_symp <- cases_adv_inc & state_df$is_symptomatic & !past_iso_delay
    to_iso  <- cases_adv_inc & state_df$is_symptomatic & past_iso_delay
    to_asymp <- cases_adv_inc & !state_df$is_symptomatic
    state_df$status[to_symp] <- "symptomatic"
    state_df$status[to_iso] <- "isolated"
    state_df$status[to_asymp] <- "asymptomatic"
    # Update rec_df
    # New symptomatic cases
    ids_to_symp <- state_df$case_id[to_symp]
    rec_df$t.symp[rec_df$case_id %in% ids_to_symp] <- sim.status$t
    rec_df$s.status[rec_df$case_id %in% ids_to_symp] <- "symptomatic"
    # New isolated cases
    ids_to_iso <- state_df$case_id[to_iso]
    rec_df$t.symp[rec_df$case_id %in% ids_to_iso] <- sim.status$t
    rec_df$t.iso[rec_df$case_id %in% ids_to_iso] <- sim.status$t
    rec_df$s.status[rec_df$case_id %in% ids_to_iso] <- "isolated"
    # New asymptomatic cases
    ids_to_asymp <- state_df$case_id[to_asymp]
    rec_df$s.status[rec_df$case_id %in% ids_to_asymp] <- "asymptomatic"
  }

  # Advance current infections from symptomatic to isolated after delay
  n_cases_adv_symp <- sum(cases_adv_symp)
  if (n_cases_adv_symp > 0){
    # Update state_df
    state_df$status[cases_adv_symp] <- "isolated"
    # Update rec_df
    ids_to_iso <- state_df$case_id[cases_adv_symp]
    rec_df$t.iso[rec_df$case_id %in% ids_to_iso] <- sim.status$t
    rec_df$s.status[rec_df$case_id %in% ids_to_iso] <- "isolated"
  }

  # Remove isolated and inactive cases from state_df
  trunc_state_df<-state_df[(state_df$status != 'isolated') & (state_df$status != 'inactive'),]

  # Add secondary infections to state_df and rec_df
  if (n_sec_infections > 0){
    out_state_df <- rbind(trunc_state_df, sec_cases_state_df)
    out_rec_df <- rbind(rec_df, sec_cases_rec_df)
  } else{
    out_state_df <- trunc_state_df
    out_rec_df <- rec_df
  }

  # Add imported infections to output state_df and rec_df
  if (n_imp_infections > 0){
    out_state_df <- rbind(out_state_df, imp_cases_state_df)
    out_rec_df <- rbind(out_rec_df, imp_cases_rec_df)
  }

  # Return updated inputs
  return(list(status=sim.status, state=out_state_df, record=out_rec_df, new_sec_cases=n_sec_infections, new_imp_cases=n_imp_infections))
}
