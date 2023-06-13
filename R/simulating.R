#' Simulating data
#'
#' Simulating data from the "threshold modulation model" (TMM) or the "parallel inhibition model" (PIM) of motor inhibition during motor imagery.
#'
#' @param nsims Numeric, number of studies to be simulated.
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param true_pars Numeric, vector of "true" parameter values.
#' @param action_mode Character, whether to simulate executed or imagined trials.
#' @param model_version Character, version of the model ("TMM" or "PIM").
#'
#' @return A dataframe containing observation (i.e., RTs and MTs).
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

simulating <- function (
        nsims = 100,
        nsamples = 2000,
        true_pars = NULL,
        action_mode = c("executed", "imagined"),
        model_version = c("TMM", "PIM")
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )
    stopifnot("true_pars must be a numeric..." = is.numeric(true_pars) )

    # testing whether only 4 pars have been specified
    stopifnot("true_pars must be a numeric of length 4..." = length(true_pars) == 4)

    # action_mode should be one of above
    action_mode <- match.arg(action_mode)

    # model_version should be one of above
    model_version <- match.arg(model_version)

    if (model_version == "TMM") {

        results <- model(
            nsims = nsims,
            nsamples = nsamples,
            exec_threshold = true_pars[4] * true_pars[1],
            imag_threshold = 0.5 * true_pars[4] * true_pars[1],
            amplitude_activ = true_pars[1],
            peak_time_activ = log(true_pars[2]),
            curvature_activ = true_pars[3],
            model_version = "TMM",
            full_output = FALSE
            ) %>%
            dplyr::mutate(action_mode = action_mode) %>%
            # keeping only the relevant columns
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
            model_version = "PIM",
            full_output = FALSE
            ) %>%
            dplyr::mutate(action_mode = action_mode) %>%
            # keeping only the relevant columns
            dplyr::select(
                .data$sim,
                reaction_time = paste0("onset_", substr(unique(.$action_mode), 1, 4) ),
                movement_time = paste0("mt_", substr(unique(.$action_mode), 1, 4) ),
                action_mode
                ) %>%
            dplyr::distinct() %>%
            dplyr::select(-.data$sim)

    }

    return (results)

}
