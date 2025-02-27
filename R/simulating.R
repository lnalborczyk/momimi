#' Simulating data
#'
#' Simulating data from the "threshold modulation model" (TMM).
#'
#' @param nsims Numeric, number of studies to be simulated.
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param true_pars Numeric, vector of "true" parameter values.
#' @param action_mode Character, whether to simulate executed or imagined trials.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param bw_noise Numeric, amount of between-trial noise.
#'
#' @return A dataframe containing observation (i.e., RTs and MTs).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # plausible "true" parameter values in the TMM
#' # parameters are the motor execution threshold, the peak time,
#' # curvature, and between-trial variability
#' true_pars <- c(1.2, 0.5, 0.5, 0.06)
#'
#' # simulating data using these parameter values
#' simulated_data <- simulating(
#'     nsims = 200,
#'     true_pars = true_pars,
#'     action_mode = "imagined"
#'     )
#'
#' # plotting the simulated data
#' plot(simulated_data)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

simulating <- function (
        nsims = 100,
        nsamples = 3000,
        true_pars = NULL,
        action_mode = c("executed", "imagined"),
        uncertainty = c("par_level", "func_level", "diffusion"),
        bw_noise = NULL
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )
    stopifnot("true_pars must be a numeric..." = is.numeric(true_pars) )

    # testing whether only 4 pars have been specified
    stopifnot("true_pars must be a numeric of 4..." = length(true_pars) == 4)

    # action_mode should be one of above
    action_mode <- match.arg(action_mode)

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    # simulating the data
    results <- model(
        nsims = nsims,
        nsamples = nsamples,
        exec_threshold = true_pars[1],
        imag_threshold = 0.5 * true_pars[1],
        peak_time = log(true_pars[2]),
        curvature = true_pars[3],
        bw_noise = true_pars[4],
        uncertainty = uncertainty,
        full_output = FALSE
        ) %>%
        dplyr::mutate(action_mode = action_mode) %>%
        dplyr::select(
            .data$sim,
            reaction_time = paste0("onset_", substr(unique(.$action_mode), 1, 4) ),
            movement_time = paste0("mt_", substr(unique(.$action_mode), 1, 4) ),
            action_mode
            ) %>%
        dplyr::distinct() %>%
        dplyr::select(-.data$sim)

    # setting the class of the resulting object
    class(results) <- c("momimi_sim", "data.frame")

    # returning it
    return (results)

}

#' @export

plot.momimi_sim <- function (x, ...) {

    optimal_bw <- x %>%
        tidyr::pivot_longer(cols = c(.data$reaction_time, .data$movement_time) ) %>%
        dplyr::pull(.data$value) %>%
        stats::bw.nrd()

    x %>%
        tidyr::pivot_longer(cols = c(.data$reaction_time, .data$movement_time) ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = .data$value,
                group = .data$name,
                colour = .data$name,
                fill = .data$name
                )
            ) +
        ggplot2::geom_histogram(
            position = "identity",
            alpha = 0.6,
            show.legend = FALSE,
            colour = "white",
            binwidth = optimal_bw
            ) +
        ggplot2::geom_label(
            data = . %>% dplyr::summarise(rt_mt_median = median(.data$value), .by = .data$name),
            ggplot2::aes(
                x = .data$rt_mt_median, y = 0,
                label = round(.data$rt_mt_median, 3),
                fill = .data$name
                ),
            position = ggplot2::position_nudge(y = 0.01),
            colour = "white",
            size = 4,
            inherit.aes = FALSE,
            show.legend = FALSE
            ) +
        ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
        ggplot2::labs(
            x = "Reaction/Movement time (s)",
            y = "Number of trials"
            )

}
