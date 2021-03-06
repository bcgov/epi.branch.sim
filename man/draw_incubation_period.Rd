% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_distributions.R
\name{draw_incubation_period}
\alias{draw_incubation_period}
\title{Draw incubation periods for new cases from distribution
defined in \code{sim_params$incub_params}}
\usage{
draw_incubation_period(n, sim_params)
}
\arguments{
\item{n}{Number of cases required}

\item{sim_params}{\code{sim_params} object (a list) containing simulation parameters, where
all of the information needed to describe the distribution is found within
\code{sim_params$incub_params} Current distributions possible are:
\itemize{
\item \code{lognormal}: Delay is drawn from lognormal distribution with
attributes "meanlog" and "sdlog" given in \code{sim_params$incub_params}.
for traced cases. Untraced cases have a further delay of "delay" days.
\item \code{weibull}: Delay is drawn from a Weibull distribution with
attributes "shape" and "scale" given in \code{sim_params$incub_params}.
Traced cases have their delay set to zero days.
}}
}
\value{
A vector of length n for case incubation period (double)
}
\description{
Draw incubation periods for new cases from distribution
defined in \code{sim_params$incub_params}
}
