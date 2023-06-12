#' Activation function
#'
#' Defining the activation/inhibition rescaled lognormal function.
#'
#' @param time Numeric, time within a trial (in seconds).
#' @param amplitude Numeric, amplitude of the activation function.
#' @param peak_time Numeric, peak time of the activation function.
#' @param curvature Numeric, curvature of the activation function.
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
        amplitude = 1.5, peak_time = 0, curvature = 0.4
        ) {

    # adding some variability in the other parameters
    # variability is currently fixed but could also be estimated
    amplitude_sim <- stats::rnorm(n = 1, mean = amplitude, sd = 0.01)
    peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = 0.01)
    curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = 0.01)

    # computing the activation/inhibition value
    activ_inhib <- amplitude_sim * exp(-(log(time) - peak_time_sim)^2 / (2 * curvature_sim^2) )

    # returning it
    return (activ_inhib)

}
