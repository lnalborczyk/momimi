% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulating.R
\name{simulating}
\alias{simulating}
\title{Simulating data}
\usage{
simulating(
  nsims = 100,
  nsamples = 5000,
  true_pars = NULL,
  action_mode = c("executed", "imagined"),
  model_version = c("TMM3", "TMM4"),
  uncertainty = c("par_level", "func_level", "diffusion"),
  bw_noise = NULL
)
}
\arguments{
\item{nsims}{Numeric, number of studies to be simulated.}

\item{nsamples}{Numeric, number of samples (time steps) within a trial.}

\item{true_pars}{Numeric, vector of "true" parameter values.}

\item{action_mode}{Character, whether to simulate executed or imagined trials.}

\item{model_version}{Character, threshold modulation model ("TMM3" or "TMM4").}

\item{uncertainty}{Numeric, indicates how noise is introduced in the system.}

\item{bw_noise}{Numeric, amount of between-trial noise.}
}
\value{
A dataframe containing observation (i.e., RTs and MTs).
}
\description{
Simulating data from the "threshold modulation model" (TMM).
}
\examples{
\dontrun{
# plausible "true" parameter values in the TMM
true_pars <- c(1.1, 0.5, 0.3, 1.25)

# simulating data using these parameter values
simulated_data <- simulating(
    nsims = 200,
    nsamples = 2000,
    true_pars = true_pars,
    action_mode = "imagined",
    model_version = "TMM4"
    )
}

}
\author{
Ladislas Nalborczyk \email{ladislas.nalborczyk@gmail.com}.
}
