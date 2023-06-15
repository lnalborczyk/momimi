#' Activation function
#'
#' Defining the activation/inhibition rescaled lognormal function.
#'
#' @param time Numeric, time within a trial (in seconds).
#' @param amplitude Numeric, amplitude of the activation function.
#' @param peak_time Numeric, peak time of the activation function.
#' @param curvature Numeric, curvature of the activation function.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param diffusion_coef Numeric, diffusion coefficient.
#' @param time_step Numeric, time step used to numerical approximation.
#'
#' @return A dataframe
#'
#' @importFrom stats rnorm
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # plotting the activation/inhibition function
#' curve(
#'     expr = activation,
#'     from = 0, 5,
#'     main = "Activation curve",
#'     xlab = "Time (in seconds)",
#'     ylab = "Activation (arbitrary units)"
#'     )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

activation <- function (
        time = 0,
        amplitude = 1.5, peak_time = 0, curvature = 0.4,
        uncertainty = c("par_specific", "overall", "brownian"),
        diffusion_coef = 0.001, time_step = 0.001
        ) {

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    if (uncertainty == "par_specific") {

        # adding some variability in the other parameters
        # variability is currently fixed but could also be estimated
        amplitude_sim <- stats::rnorm(n = 1, mean = amplitude, sd = 0.01)
        peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = 0.01)
        curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = 0.01)

        activ_inhib <- amplitude_sim *
            exp(-(log(time) - peak_time_sim)^2 / (2 * curvature_sim^2) )

    } else if (uncertainty == "overall") {

        activ_inhib <- pmax(amplitude *
            exp(-(log(time) - peak_time)^2 / (2 * curvature^2) ) +
            stats::rnorm(n = 1, mean = 0, sd = 0.01), 0)

    } else if (uncertainty == "brownian") {

        # derivative of activation
        d_activation <- function (t, A, mu, sigma) {

            x <- -((log(t) - mu) / (t * sigma^2) ) *
                (A * exp(-((log(t) - mu)^2) / (2 * sigma^2) ) )

            return (x)

        }

        # diffusion process
        activ_inhib <- time_step * d_activation(
            t = time,
            A = amplitude, mu = peak_time, sigma = curvature
            ) + diffusion_coef * sqrt(time_step) * stats::rnorm(n = 1, mean = 0, sd = 1)

        # cumulative sum
        activ_inhib <- pmax(cumsum(tidyr::replace_na(data = activ_inhib, replace = 0) ), 0)

    }

    # returning it
    return (activ_inhib)

}
