% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss.R
\name{loss}
\alias{loss}
\title{Loss function}
\usage{
loss(
  par,
  data,
  nsims = NULL,
  nsamples = 3000,
  exec_threshold = 1,
  imag_threshold = 0.5,
  uncertainty = c("par_level", "func_level", "diffusion"),
  time_step = 0.001,
  diffusion_coef = 0.001,
  error_function = c("g2", "rmse", "sse", "wsse", "ks")
)
}
\arguments{
\item{par}{Numeric, vector of parameter values.}

\item{data}{Dataframe, data to be used for computing the error.}

\item{nsims}{Numeric, number of simulations (observations/trials).}

\item{nsamples}{Numeric, number of samples (time steps) within a trial.}

\item{exec_threshold}{motor execution threshold.}

\item{imag_threshold}{motor imagery threshold.}

\item{uncertainty}{Numeric, indicates how noise is introduced in the system.}

\item{time_step}{Numeric, time step used to numerical approximation.}

\item{diffusion_coef}{Numeric, diffusion coefficient.}

\item{error_function}{Character, loss function to be used when fitting the model.}
}
\value{
The loss/error for a given set of observations.
}
\description{
Provides the loss function for the model.
}
\author{
Ladislas Nalborczyk \email{ladislas.nalborczyk@gmail.com}.
}
