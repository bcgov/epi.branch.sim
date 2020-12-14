epi.branch.sim
==============

[![img](https://img.shields.io/badge/Lifecycle-Stable-97ca00)](https://github.com/bcgov/repomountie/blob/master/doc/lifecycle-badges.md)

Simulates an epidemic outbreak with a stochastic branching process from
a number of initial seed cases.

### Features

Every simulated case has disease and intervention parameters
stochastically drawn from user-provided distributions. These parameters
determine which secondary cases occur, which in turn create the next
generation of cases. This simulation implements modelling of
interventions such as: \* manual/conventional contact tracing, \*
app-based contact tracing usage and compliance, \* reduction of contacts
through physical distancing, \* self-isolation of cases after certain
conditions are met, and \* quarantine/isolation of new arrivals.

This simulation also includes the option of adding new cases via
importation, which may have different contact tracing and self-isolation
conditions. The simulation itself outputs a dataframe that summarizes
every case generated during the course of the simulation. This can be
used to generate metrics or other statistics across a large number of
runs. Some example driver functions to do this are provided in the
package and as vignettes.

This package is inspired by [Hellewell et al. “Feasibility of
controlling COVID-19 outbreaks by isolation of cases and contacts.” *The
Lancet Global Health* 2020; 8: E488–E496. DOI:
10.1016/S2214-109X(20)30074-7](https://doi.org/10.1016/S2214-109X(20)30074-7).

### Installation

To install this package along with the vignettes, use

``` r
remotes::install_github("bcgov/epi.branch.sim", ref="main", build_vignettes=TRUE)
```

### Usage

A very simple example of setting up and running a scenario for 30 days.

``` r
# Set up simulation objects
sim_params <- epi.branch.sim::initialize_sim_params(
  R0=3.0, infect_dur=999, vary_trace=FALSE, p_trace=0.8, 
  p_trace_app=0, p_trace_app_comp=0, p_symp=0.9, dt=1,
  incub_params=list(dist='weibull', shape=2.322737, scale=6.492272),
  generation_int_params=list(dist='skew_norm', omega=2, alpha=1.95),
  iso_delay_params=list(dist='Hellewell', shape=1.651524, scale=4.287786),
  sec_infect_params=list(type='Hellewell', disp=0.16),
  import_params=list(type="None"), 
  phys_dist_params=list(
    pd_pop_frac = 0, 
    pd_contact_rate1 = 1.0,
    pd_contact_rate2 = 1.0,
    pd_change_t = 0)
)
start_time <- 0
n_initial <- 20
sim_status <- epi.branch.sim::initialize_sim_status(start_time,n_initial) 
state_df   <- epi.branch.sim::create_state_df(n_initial,sim_params, sim_status, initialize=TRUE)
record_df  <- epi.branch.sim::create_record_df(state_df, sim_status, initialize=TRUE)

# Run for 30 steps
for (t in 1:30){
  out <- epi.branch.sim::step_simulation(sim_status, state_df, record_df, sim_params)
  sim_status <- out$status # update sim_status
  state_df <- out$state    # update state_df
  record_df <- out$record  # update record_df
}

# Output products are the state_df and record_df data frames
```

Detailed documentation on the simulation algorithm, input parameters,
output data frames and examples are provided as vignettes. Once the
package is installed with the `build_vignette=TRUE` option as above,
these vignettes can be viewed with:

``` r
browseVignettes('epi.branch.sim')
```

For help on all documented functions and objects, use:

``` r
help(package='epi.branch.sim')
```

### Project Status

This package is stable. This is a public release of a tool regularly
used within our team in order to share the tool more widely. At release,
the package is functional, all of the intended features are working and
documented. More features may be added based on usage and demand, but there
are no current plans to continue development at this time.

### Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an
[issue](https://github.com/bcgov/epi.branch.sim/issues/).

### How to Contribute

If you would like to contribute to the package, please see our
[CONTRIBUTING](CONTRIBUTING.md) guidelines.

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.

### License

    Copyright 2020 Province of British Columbia

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software 
    distributed under the License is distributed on an "AS IS" BASIS, 
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and 
    limitations under the License.

------------------------------------------------------------------------

*This project was created using the
[bcgovr](https://github.com/bcgov/bcgovr) package.*
