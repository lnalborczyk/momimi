#' Activation function
#'
#' Defining the activation/inhibition rescaled lognormal function.
#'
#' @param time Numeric, time within a trial (in seconds).
#' @param peak_time Numeric, peak time of the activation function.
#' @param curvature Numeric, curvature of the activation function.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param bw_noise Numeric, amount of between-trial noise.
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
#'     ylab = "Activation (a.u.)"
#'     )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

activation <- function (
        time = 0,
        peak_time = 0.5, curvature = 0.4,
        uncertainty = c("par_level", "func_level", "diffusion"),
        bw_noise = 0.1, diffusion_coef = 0.001, time_step = 0.001
        ) {

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    if (uncertainty == "par_level") {

        peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = bw_noise)
        # curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = bw_noise)

        # activ_inhib <- exp(-(log(time) - peak_time_sim)^2 / (2 * curvature_sim^2) )
        activ_inhib <- exp(-(log(time) - peak_time_sim)^2 / (2 * curvature^2) )

    } else if (uncertainty == "func_level") {

        activ_inhib <- pmax(
            exp(-(log(time) - peak_time)^2 / (2 * curvature^2) ) +
            stats::rnorm(n = 1, mean = 0, sd = bw_noise), 0
            )

    } else if (uncertainty == "diffusion") {

        # derivative of activation function
        d_activation <- function (t, mu, sigma) {

            # x <- -((log(t) - mu) / (t * sigma^2) ) *
            #     (A * exp(-((log(t) - mu)^2) / (2 * sigma^2) ) )

            x <- -(exp(-(log(t) - mu)^2 / (2 * sigma^2) ) *
                       (2 * (1 / t * (log(t) - mu) ) / (2 * sigma^2) ) )

            return (x)

        }

        # diffusion process
        activ_inhib <- time_step * d_activation(
            t = time,
            mu = peak_time, sigma = curvature
            ) + diffusion_coef * sqrt(time_step) * stats::rnorm(n = 1, mean = 0, sd = 1)

        # replacing the initial NA by 0
        activ_inhib <- tidyr::replace_na(data = activ_inhib, replace = 0)

        # cumulative sum (absorbing boundary)
        activ_inhib <- purrr::accumulate(.x = activ_inhib, .f = ~ ifelse(.x + .y < 0, 0, .x + .y) )

    }

    # returning it
    return (activ_inhib)

}
