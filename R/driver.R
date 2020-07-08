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

# driver.R
# Contains the functions that will create simulation objects, run a number of simulations
# corresponding to a given scenario, and then tabulate metrics and summary statistics to
# be used in analysis. The functions provided here are meant to be examples and/or
# convienence functions. Users of this package may want to write their own driver
# functions to suit their own analysis needs.

#' Run scenarios as in Hellewell et al. (2020) paper
#'
#' For a given tracing fraction and R0 (a "scenario"), count the fraction of simulations
#' that are controlled. Designed to reproduce results from the source paper for this model: [Hellewell et al. "Feasibility of controlling COVID-19 outbreaks by isolation of cases and contacts." *The Lancet Global Health* 2020; 8: E488-E496. DOI: 10.1016/S2214-109X(20)30074-7](https://doi.org/10.1016/S2214-109X(20)30074-7)
#'
#' The intention of this driver function is to reproduce Figure 3A from Hellewell et al. which
#' shows the fraction of outbreaks controlled as a function of R0 and tracing fraction. These
#' parameters are the positional arguments to this function. The remaining named arguments allow
#' for consideration of other scenarios examined by Hellewell et al. The named arguments' default
#' values are Hellewell et al.'s benchmark values but can be changed to reproduce results from
#' Figure 4 of the source paper.
#'
#' The secondary intention of this driver function is to provide a simple example of how to
#' call the exported functions of this module in order to set up simulations objects, run
#' a number of simulations and return summary statistics/metrics for further analysis.
#'
#' @param nsims     Number of simulations to run in this scenario
#' @param p_trace   Fraction of new infections that are traced
#' @param R0        Basic reproduction number for disease
#' @param infections_before_symptoms  One of "low", "medium" or "high" corresponding to one
#'                                    of three generation interval distributions considered by
#'                                    Hellewell et al. These distributions result in <1%,
#'                                    15%, and 30% of generation interval lengths occuring prior
#'                                    to symptom onset. "Medium" (15%) is the default value
#'                                    and the benchmark used in Hellewell et al.
#' @param iso_delay_length  One of "short" or "long" corresponding to two distributions used
#'                          by Hellewell for the delay from symptom onset to isolation for
#'                          untraced symptomatic cases. "Short" is the default and
#'                          benchmark value.
#' @param p_symp    Fraction of cases that are symptomatic. 1 is the default and benchmark value.
#' @param initial_cases   Number of cases at the start. 20 is the default and benchmark value
#' @param max_cases       Number of total cases before simulation ends with an "epidemic"
#'                        outcome. Defaults to 5000, used by Hellewell et al.
#'
#' @return A named vector with the following entries describing the scenario and
#' the outcome (percentage of simulations that went extinct and absolute number
#' of each type of outcome). Namely the entries are: nsims, R0, p_trace,
#' percentage, nsims_extinct, nsims_epidemic, nsims_ongoing.
#' @export
run_scenario_Hellewell <- function(nsims,p_trace,R0,
                                   infections_before_symptoms='medium',
                                   iso_delay_length='short',
                                   p_symp = 1.0,
                                   initial_cases = 20,
                                   max_cases = 5000){
  col_names <- c("nsims", "R0", "p_trace","percent_controlled", "n_extinct", "n_epidemic", "n_ongoing")
  cat(sprintf('Running %.0f simulations of scenario: R0 = %.2f, p_trace = %.2f\n',nsims,R0,p_trace))
  # Counters for simulation outcomes
  nsims_extinct<-0 # all cases went extinct
  nsims_epidemic<-0  # more than [max_cases] total cases
  nsims_ongoing<-0   # Less than [max_cases] total cases but still new infections after 12 weeks
  for (j in 1:nsims){
    # Set up incubation period length distribution as per Hellewell et al.
    incub_params <- list(dist='weibull', shape=2.322737, scale=6.492272)

    # Set up generation intervals to one of three levels used by Hellewell et al. (called serial intervals there)
    if (infections_before_symptoms == 'low'){ # <1% infections before symptom onset
      generation_int_params <- list(dist='skew_norm', omega=2, alpha=30)
    } else if (infections_before_symptoms == 'medium'){ # 15% infections before symptom onset
      generation_int_params <- list(dist='skew_norm', omega=2, alpha=1.95) # DEFAULT BENCHMARK
    } else if (infections_before_symptoms == 'high'){ # 30% infections before symptoms
      generation_int_params <- list(dist='skew_norm', omega=2, alpha=0.7)
    } else {
      stop('Invalid value for argument "infections_before_symptoms".')
    }

    # Set up isolation delay parameters to one of two levels used by Hellewell et al.
    if (iso_delay_length == 'short'){
      iso_delay_params <- list(dist='Hellewell', shape=1.651524, scale=4.287786) # BENCHMARK
    } else if (iso_delay_length == 'long'){
      iso_delay_params <- list(dist='Hellewell', shape=2.305172, scale=9.483875)
    } else{
      stop('Invalid value for argument "iso_delay_length".')
    }

    # Set up dispersion parameter in number of secondary infections
    sec_infect_params <- list(type='Hellewell', disp=0.16) # From JH2020 "overdispersion in R0"

    # Turn off / disable features not used by Hellewell's work
    infect_dur <- 999 # Hellewell et al. doesn't deactivate cases so we won't either
    do_variable_trace<-FALSE
    p_trace_app<-0
    p_trace_app_comp<-0
    dt <- 1 # Advance one day at a time
    import_params <- list(type='None')
    phys_dist_params <- list(pd_pop_frac=0,pd_contact_rate1=1,pd_contact_rate2=1,pd_change_t=0)


    # Initialize the required simulation objects
    sim_params <- initialize_sim_params(
      R0, infect_dur, do_variable_trace, p_trace,
      p_trace_app, p_trace_app_comp, p_symp, dt,
      incub_params, generation_int_params,
      iso_delay_params, sec_infect_params,
      import_params, phys_dist_params
    )
    sim_status <- initialize_sim_status(start_time=0,initial_cases)
    state_df   <- create_state_df(initial_cases,sim_params,sim_status,initialize=TRUE)
    record_df  <- create_record_df(state_df, sim_status, initialize=TRUE)

    timemax=ceiling(12*7/dt) # 12 weeks after initial cases (Fig 3 caption)
    n_new_cases=rep(0,timemax) # number of new cases per step

    early_exit <- FALSE
    # Run simulation
    for (ii in 1:timemax){
      out <- step_simulation(sim_status, state_df, record_df, sim_params)
      sim_status <- out$status
      state_df <- out$state
      record_df <- out$record
      n_new_cases[ii] <-out$new_sec_cases
      if (nrow(state_df)==0){
        # stop running if no more contagious cases
        nsims_extinct <- nsims_extinct + 1
        early_exit <- TRUE
        break
      } else if(nrow(record_df)>max_cases){
        # also stop running if more than [max_cases] cases (not controlled)
        nsims_epidemic <- nsims_epidemic + 1
        early_exit <- TRUE
        break
      }
    }
    # For runs that finished running, there must still be cases left, so they are "ongoing"
    if (!early_exit){
      if (nrow(state_df)>0){
        nsims_ongoing <- nsims_ongoing+1
      }
    }
  }

  # Create output
  percentage = nsims_extinct/nsims * 100
  results <- c(nsims, R0, p_trace, percentage,nsims_extinct,nsims_epidemic,nsims_ongoing)
  names(results)<- col_names

  return(results)
}

# ABC_start = sim day number to begin ABC matching
# ABC_end = sim day number to end ABC matching
# ABC_match_type = "reported" if using real reports or "adjusted" if adjusting for delay
# ABC_match_last_n = match to last n rows of BC case data file
# ABC_region = "BC" or "VIHA"
run_scenarios <- function(scenario_params, outdir='.',
                          do_ABC=FALSE, ABC_start=7, ABC_end=13,
                          ABC_match_type='reported', ABC_match_last_n=7,
                          ABC_region='BC'){
  scn_id <- scenario_params$scn_id
  nsims <- scenario_params$nsims
  dt <- scenario_params$dt
  tmax_days <- scenario_params$tmax_days
  max_active_cases <- scenario_params$max_active_cases
  save_scn <- scenario_params$save_scn
  import_model <- scenario_params$import_model
  do_variable_trace <- scenario_params$do_variable_trace
  contact_rates_str <- scenario_params$contact_rates
  contact_rates <- eval(parse(text=contact_rates_str))
  pd_p_group_str <- scenario_params$pd_p_group
  pd_p_group <- eval(parse(text=pd_p_group_str))
  pd_delay <- scenario_params$pd_delay
  R0 <- scenario_params$R0
  p_trace <- scenario_params$p_trace
  p_symp <- scenario_params$p_symp
  n_initial <- scenario_params$n_initial

  results=NULL
  for (j in 1:nsims){
    infect_dur <- 999 # Hellewell et al. doesn't deactivate cases so we won't either
    incub_params <- list(dist='lognormal',meanlog=1.57, sdlog=0.65)
    generation_int_params <- list(dist='gamma', shape=2.29, rate=0.36)
    iso_delay_params <- list(dist='uniform_delay', min=2, max=3, delay=2) # 2-3 days if traced, 4-5 if not
    sec_infect_params <- list(type='Hellewell', disp=0.58) # From JH2020 but with larger disp.
    if (import_model=='None'){
      import_params <- 'None'
    }
    else{
      import_params <- get_import_params(file.path('data',paste0(import_model,'.csv')))
    }
    phys_dist_params <- list(contact_rates=contact_rates, p_group=pd_p_group, delay=pd_delay)
    sim_params <- initialize_sim_params(R0, infect_dur, do_variable_trace, p_trace,
                                        p_symp, dt,
                                        incub_params, generation_int_params,
                                        iso_delay_params, sec_infect_params,
                                        import_params, phys_dist_params)
    sim_status <- initialize_sim_status(0,n_initial)
    state_df   <- create_state_df(n_initial, sim_params, sim_status, initialize=TRUE)
    record_df  <- create_record_df(state_df, sim_status, initialize=TRUE)

    timemax=ceiling(tmax_days/dt)
    # early_exit <- FALSE
    # Set up vectors to track key metrics per time step
    n_total=rep(NA,timemax) # number of total cases
    n_active=rep(NA,timemax) # number of active cases
    n_incub=rep(NA,timemax) # number of cases incubating
    n_symp=rep(NA,timemax) # number of cases with symptoms
    n_asymp=rep(NA,timemax) # number of cases not showing symptoms but past incubation
    n_iso=rep(NA,timemax) # number of cases isolated
    n_new_sec=rep(NA,timemax) # number of new secondary cases
    n_new_imp=rep(NA,timemax) # number of new imported cases
    R0_eff=rep(NA,timemax)
    reject_run<-FALSE
    for (ii in 1:timemax){
      # Take a step forward
      out <- step_simulation(sim_status, state_df, record_df, sim_params)
      sim_status <- out$status
      state_df <- out$state
      record_df <- out$record
      # Track key metrics
      n_total[ii] <-nrow(record_df) # all cases ever
      n_active[ii] <-nrow(state_df) # incub + sympt + asympt
      n_incub[ii] <- sum(state_df$status=='incubation')
      n_symp[ii] <- sum(state_df$status=='symptomatic')
      n_asymp[ii] <- sum(state_df$status=='asymptomatic')
      n_iso[ii] <- sum(record_df$s_status=='isolated') # get this from record_df!
      n_new_sec[ii] <-out$new_sec_cases # new sec cases in last dt
      n_new_imp[ii] <- out$new_imp_cases # new imp cases in last dt
      R0_eff[ii] <- mean(record_df$n_sec_infects)
      # If doing ABC, check for match on last day of ABC matching
      if (do_ABC & ii==ABC_end/dt){
        reported_counts <- get_BC_counts(region=ABC_region,last_n=ABC_match_last_n,match_type=ABC_match_type)
        under_report_fac <- scenario_params$under_report_fac # e.g. 4x
        corrected_counts <- reported_counts * under_report_fac
        if (under_report_fac > 1){
          sim_cases <- n_total[seq(ABC_start/dt,ABC_end/dt,1/dt)] # Match with total case counts
        } else{
          sim_cases <- n_iso[seq(ABC_start/dt,ABC_end/dt,1/dt)] # Match with diagnosed cases
        }
        is_under <- any(sim_cases < corrected_counts*0.75) # boolean: any row below 75%?
        is_over <- any(sim_cases > corrected_counts*1.25) # boolean: any row above 125%
        if (is_under | is_over){ # not a match
          reject_run<-TRUE # flag to not save this simulation run into the output df
          break # end sim early
        }
      }
    }
    # Collapse metrics into daily counts (right now assuming dt=1, need a way to sum it up)
    day=0:(timemax*dt)
    bin_every = 1/dt
    nd_total <- c(n_initial,n_total[seq(0,timemax,bin_every)])
    nd_active <- c(n_initial,n_active[seq(0,timemax,bin_every)])
    nd_incub <- c(n_initial,n_incub[seq(0,timemax,bin_every)])
    nd_symp <- c(0,n_symp[seq(0,timemax,bin_every)])
    nd_asymp <- c(0,n_asymp[seq(0,timemax,bin_every)])
    nd_iso <- c(0,n_iso[seq(0,timemax,bin_every)])
    nd_new_S <- c(0,colSums(matrix(n_new_sec,nrow=bin_every)))
    nd_new_I <- c(0,colSums(matrix(n_new_imp,nrow=bin_every)))
    nd_R0eff <- c(R0,R0_eff[seq(0,timemax,bin_every)])
    # Create output dataframe for this run
    results_sim <- data.frame(
      # scenario parameters
      "R0" = R0,
      "p_trace" = p_trace,
      "p_symp" = p_symp,
      "initial_n" = n_initial,
      "sim_id" = sprintf('%.0f.%.0f',scn_id,j), # scenario #, sim #
      # results
      "day" = day,
      "n_total" = nd_total,
      "n_active" = nd_active,
      "n_incub" = nd_incub,
      "n_symp" = nd_symp,
      "n_asymp" = nd_asymp,
      "n_iso" = nd_iso,
      "n_new_S" = nd_new_S,
      "n_new_I" = nd_new_I,
      "R0eff" = nd_R0eff
      )
    # Append to dataframe
    if (!reject_run){
      results <- rbind(results,results_sim)
    }
    # Save a temporary dataframe after each simulation
    if (save_scn){
      s_fname=sprintf('df_scn%04.0f_temp.rds',scn_id)
      saveRDS(results,file=file.path(outdir,'temp',s_fname))
    }
  }
  # Save a final dataframe at the end of the scenario
  if (save_scn){
    s_fname=sprintf('df_scn%04.0f.rds',scn_id)
    saveRDS(results,file=file.path(outdir,s_fname))
  }

  return(results)
}
