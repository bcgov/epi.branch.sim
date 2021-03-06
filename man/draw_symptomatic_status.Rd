% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_distributions.R
\name{draw_symptomatic_status}
\alias{draw_symptomatic_status}
\title{Determine whether a new case is symptomatic from a uniform distribution}
\usage{
draw_symptomatic_status(n, sim_params)
}
\arguments{
\item{n}{Number of cases required}

\item{sim_params}{\code{sim_params} object (a list) containing simulation parameters, here,
the \code{sim_params$p_sym} value is used as the probability of a case
being symptomatic.}
}
\value{
A boolean vector of length n for whether a case is symptomatic
}
\description{
Determine whether a new case is symptomatic from a uniform distribution
}
