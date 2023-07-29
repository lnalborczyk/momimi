#' Simulating data
#'
#' Simulating data from the "threshold modulation model" (TMM) or the "parallel inhibition model" (PIM) of motor inhibition during motor imagery.
#'
#' @param nsims Numeric, number of studies to be simulated.
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param true_pars Numeric, vector of "true" parameter values.
#' @param action_mode Character, whether to simulate executed or imagined trials.
#' @param model_version Character, threshold modulation model ("TMM3" or "TMM4") or parallel inhibition model ("PIM").
#'
#' @return A dataframe containing observation (i.e., RTs and MTs).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # plausible "true" parameter values in the TMM
#' true_pars <- c(1.1, 0.5, 0.3, 1.25)
#'
#' # simulating data using these parameter values
#' simulated_data <- simulating(
#'     nsims = 200,
#'     nsamples = 2000,
#'     true_pars = true_pars,
#'     action_mode = "imagined",
#'     model_version = "TMM4"
#'     )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

simulating <- function (
        nsims = 100,
        nsamples = 2000,
        true_pars = NULL,
        action_mode = c("executed", "imagined"),
        model_version = c("TMM3", "TMM4", "PIM")
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )
    stopifnot("true_pars must be a numeric..." = is.numeric(true_pars) )

    # testing whether only 3 or 4 pars have been specified
    stopifnot("true_pars must be a numeric of length 3 or 4..." = length(true_pars) %in% c(3, 4) )

    # action_mode should be one of above
    action_mode <- match.arg(action_mode)

    # model_version should be one of above
    model_version <- match.arg(model_version)

    if (model_version == "TMM3") {

        results <- model(
            nsims = nsims,
            nsamples = nsamples,
            exec_threshold = true_pars[1] * 1.5,
            imag_threshold = 0.5 * true_pars[1] * 1.5,
            amplitude_activ = 1.5,
            peak_time_activ = log(true_pars[2]),
            curvature_activ = true_pars[3],
            model_version = model_version,
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

    } else if (model_version == "TMM4") {

        results <- model(
            nsims = nsims,
            nsamples = nsamples,
            exec_threshold = true_pars[4] * true_pars[1],
            imag_threshold = 0.5 * true_pars[4] * true_pars[1],
            amplitude_activ = true_pars[1],
            peak_time_activ = log(true_pars[2]),
            curvature_activ = true_pars[3],
            model_version = model_version,
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

    } else if (model_version == "PIM") {

        results <- model(
            nsims = nsims,
            nsamples = nsamples,
            exec_threshold = 1,
            imag_threshold = 0.5,
            amplitude_activ = 1.5,
            peak_time_activ = log(true_pars[2]),
            curvature_activ = true_pars[3],
            amplitude_inhib = 1.5 / true_pars[1],
            peak_time_inhib = log(true_pars[2]),
            curvature_inhib = true_pars[4] * true_pars[3],
            model_version = model_version,
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

    }

    # setting the class of the resulting object
    class(results) <- c("momimi_sim", "data.frame")

    # returning it
    return (results)

}

plot.momimi_sim <- function (x, ...) {

    x %>%
        tidyr::pivot_longer(cols = c(.data$reaction_time, .data$movement_time) ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = .data$value, group = .data$name,
                colour = .data$name, fill = .data$name
                )
            ) +
        ggplot2::geom_density(
            color = "white",
            alpha = 0.6,
            adjust = 5,
            show.legend = FALSE
            ) +
        ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
        ggplot2::labs(
            x = "Reaction/Movement time (in seconds)",
            y = "Probability density"
            )

}
