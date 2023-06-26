#' Generating initial population
#'
#' Generating plausible initial parameter values for DEoptim based on empirical constraints.
#'
#' @param nstudies Numeric, number of starting values in the LHS.
#' @param action_mode Character, action mode (executed or imagined).
#' @param par_names Character, vector of parameter names.
#' @param lower_bounds Numeric, vector of lower bounds for parameters.
#' @param upper_bounds Numeric, vector of upper bounds for parameters.
#' @param model_version Character, threshold modulation model ("TMM") or parallel inhibition model ("PIM").
#' @param rt_contraints Numeric, vector of length 2 specifying the min and max RT (in seconds).
#' @param mt_contraints Numeric, vector of length 2 specifying the min and max MT (in seconds).
#' @param verbose Boolean, whether to print progress during fitting.
#'
#' @return Returns a dataframe containing plausible (according to custom constraints) initial parameter values.
#'
#' @importFrom rlang .data
#' @importFrom stats rnorm
#' @importFrom magrittr %>%
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

generating_initialpop <- function (
        nstudies, action_mode,
        par_names, lower_bounds, upper_bounds,
        model_version = c("TMM", "PIM"),
        rt_contraints = c(0.1, 2), mt_contraints = c(0.1, 2),
        verbose = TRUE
        ) {

    # some tests for variable types
    stopifnot("nstudies must be a numeric..." = is.numeric(nstudies) )
    stopifnot("action_mode must be a character..." = is.character(action_mode) )

    # testing whether only 4 pars have been specified
    stopifnot("par_names must be a numeric of length 4..." = length(par_names) == 4)
    stopifnot("lower_bounds must be a numeric of length 4..." = length(lower_bounds) == 4)
    stopifnot("upper_bounds must be a numeric of length 4..." = length(upper_bounds) == 4)
    stopifnot("rt_contraints must be a numeric of length 2..." = length(rt_contraints) == 2)
    stopifnot("mt_contraints must be a numeric of length 2..." = length(mt_contraints) == 2)

    # model_version should be one of above
    model_version <- match.arg(model_version)

    # initialising the result_nrow variable
    result_nrow <- 0

    # while we do not have the requested number of starting values
    while (result_nrow < nstudies) {

        # initialising an empty dataframe
        lhs_pars <- data.frame(matrix(data = NA, nrow = nstudies, ncol = length(lower_bounds) ) )

        # populating it with hypercube samples
        for (i in 1:length(lower_bounds) ) {

            lhs_pars[, i] <- tgp::lhs(n = nstudies, rect = c(lower_bounds[i], upper_bounds[i]) )[, 1]

        }

        # setting columns names
        colnames(lhs_pars) <- par_names

        if (model_version == "TMM") {

            # defining the activation/inhibition rescaled lognormal function
            activation_function <- function (
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude = 1.5, peak_time = 0, curvature = 0.4
        ) {

                # adding some variability in the other parameters
                # variability is currently fixed but could also be estimated
                amplitude_sim <- rnorm(n = 1, mean = amplitude, sd = 0.01)
                peak_time_sim <- rnorm(n = 1, mean = peak_time, sd = 0.01)
                curvature_sim <- rnorm(n = 1, mean = curvature, sd = 0.01)
                exec_threshold_sim <- rnorm(n = 1, mean = exec_threshold, sd = 0.01)

                # # no variability in the motor imagery threshold
                # imag_threshold_sim <- rnorm(n = 1, mean = imag_threshold, sd = 0.01)
                imag_threshold_sim <- imag_threshold

                # computing the predicted RT and MT in imagery
                onset_offset_imag <- onset_offset(
                    amplitude_activ = amplitude_sim,
                    peak_time_activ = peak_time_sim,
                    curvature_activ = curvature_sim,
                    thresh = imag_threshold_sim,
                    model_version = model_version
                    )

                onset_imag <- min(onset_offset_imag)
                mt_imag <- max(onset_offset_imag) - min(onset_offset_imag)

                # computing the predicted RT and MT in execution
                # computing the predicted RT and MT in imagery
                onset_offset_exec <- onset_offset(
                    amplitude_activ = amplitude_sim,
                    peak_time_activ = peak_time_sim,
                    curvature_activ = curvature_sim,
                    thresh = exec_threshold_sim,
                    model_version = model_version
                    )

                onset_exec <- min(onset_offset_exec)
                mt_exec <- max(onset_offset_exec) - min(onset_offset_exec)

                # returning it
                return (data.frame(onset_imag, mt_imag, onset_exec, mt_exec) )

            }

            # computing the predicted RT and MT
            predicted_rt_mt <- lhs_pars %>%
                dplyr::rowwise() %>%
                dplyr::do(
                    suppressWarnings(
                        activation_function(
                            exec_threshold = .data$exec_threshold * .data$amplitude_activ,
                            imag_threshold = 0.5 * .data$exec_threshold * .data$amplitude_activ,
                            amplitude = .data$amplitude_activ,
                            peak_time = log(.data$peak_time_activ),
                            curvature = .data$curvature_activ
                            )
                        )
                    )

        } else if (model_version == "PIM") {

            # defining the balance function
            # basically a ratio of two rescaled lognormal functions
            balance_function <- function (
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
        amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6
        ) {

                # adding some variability in the other parameters
                # variability is currently fixed but could also be estimated
                amplitude_activ_sim <- stats::rnorm(n = 1, mean = amplitude_activ, sd = 0.01)
                peak_time_activ_sim <- stats::rnorm(n = 1, mean = peak_time_activ, sd = 0.01)
                curvature_activ_sim <- stats::rnorm(n = 1, mean = curvature_activ, sd = 0.01)

                amplitude_inhib_sim <- stats::rnorm(n = 1, mean = amplitude_inhib, sd = 0.01)
                peak_time_inhib_sim <- stats::rnorm(n = 1, mean = peak_time_inhib, sd = 0.01)
                curvature_inhib_sim <- stats::rnorm(n = 1, mean = curvature_inhib, sd = 0.01)

                # computing the predicted RT and MT in imagery
                onset_offset_imag <- onset_offset(
                    amplitude_activ = amplitude_activ_sim,
                    peak_time_activ = peak_time_activ_sim,
                    curvature_activ =  curvature_activ_sim,
                    amplitude_inhib = amplitude_inhib_sim,
                    peak_time_inhib = peak_time_inhib_sim,
                    curvature_inhib = curvature_inhib_sim,
                    thresh = imag_threshold,
                    model_version = model_version
                    )

                onset_imag <- min(onset_offset_imag)
                mt_imag <- max(onset_offset_imag) - min(onset_offset_imag)

                # computing the predicted RT and MT in execution
                onset_offset_exec <- onset_offset(
                    amplitude_activ = amplitude_activ_sim,
                    peak_time_activ = peak_time_activ_sim,
                    curvature_activ =  curvature_activ_sim,
                    amplitude_inhib = amplitude_inhib_sim,
                    peak_time_inhib = peak_time_inhib_sim,
                    curvature_inhib = curvature_inhib_sim,
                    thresh = exec_threshold,
                    model_version = model_version
                    )

                onset_exec <- min(onset_offset_exec)
                mt_exec <- max(onset_offset_exec) - min(onset_offset_exec)

                # returning it
                return (data.frame(onset_imag, mt_imag, onset_exec, mt_exec) )

            }

            # computing the predicted RT and MT
            predicted_rt_mt <- lhs_pars %>%
                dplyr::rowwise() %>%
                dplyr::do(
                    suppressWarnings(
                        balance_function(
                            amplitude_activ = 1.5,
                            peak_time_activ = log(.data$peak_time),
                            curvature_activ = .data$curvature_activ,
                            amplitude_inhib = 1.5 / .data$amplitude_ratio,
                            peak_time_inhib = log(.data$peak_time),
                            curvature_inhib = .data$curvature_inhib * .data$curvature_activ
                            )
                        )
                    )

        }

        if (action_mode == "imagined") {

            predicted_rt_mt <- predicted_rt_mt[, 1:2]

        } else if (action_mode == "executed") {

            predicted_rt_mt <- predicted_rt_mt[, 3:4]

        }

        # adding some constraints
        # to restrict parameter values to realistic RT/MT data (e.g., Evans, 2020)
        # predicted RT/MT should be valid (not a NaN)
        # predicted RT should be be between 0.2 and 1 seconds
        # predicted MT should be be between 0.2 and 2 seconds
        # balance at the end of the trial should come back below 0.25
        if (model_version == "TMM") {

            final_par_values <- dplyr::bind_cols(lhs_pars, predicted_rt_mt) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    balance_end_of_trial = .data$amplitude_activ *
                        exp(-(log(3) - .data$peak_time_activ)^2 / (2 * .data$curvature_activ^2) )
                    ) %>%
                dplyr::mutate(
                    included = dplyr::case_when(
                        any(is.na(dplyr::pick(dplyr::everything() ) ) ) ~ FALSE,
                        dplyr::pick(length(par_names) + 1) < min(rt_contraints) ~ FALSE,
                        dplyr::pick(length(par_names) + 1) > max(rt_contraints) ~ FALSE,
                        dplyr::pick(length(par_names) + 2) < min(mt_contraints) ~ FALSE,
                        dplyr::pick(length(par_names) + 2) > max(mt_contraints) ~ FALSE,
                        balance_end_of_trial > 0.25 * exec_threshold * amplitude_activ ~ FALSE,
                        .default = TRUE
                        )
                    )

        } else if (model_version == "PIM") {

            final_par_values <- dplyr::bind_cols(lhs_pars, predicted_rt_mt) %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    balance_end_of_trial = .data$amplitude_ratio *
                        exp(-(log(3) - .data$peak_time)^2 / (2 * .data$curvature_activ^2) +
                                (log(3) - .data$peak_time)^2 / (2 * .data$curvature_inhib^2) )
                    ) %>%
                dplyr::mutate(
                    included = dplyr::case_when(
                    any(is.na(dplyr::pick(dplyr::everything() ) ) ) ~ FALSE,
                    dplyr::pick(length(par_names) + 1) < min(rt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 1) > max(rt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 2) < min(mt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 2) > max(mt_contraints) ~ FALSE,
                    balance_end_of_trial > 0.25 ~ FALSE,
                    .default = TRUE
                    ) )

        }

        # reshaping the final matrix
        temp_lhs_initial_pop <- final_par_values %>%
            dplyr::filter(.data$included) %>%
            dplyr::select(1:length(par_names) ) %>%
            data.frame() %>%
            as.matrix()

        # rbinding results
        if (!exists("lhs_initial_pop") ) {

            lhs_initial_pop <- temp_lhs_initial_pop

        } else {

            lhs_initial_pop <- rbind(lhs_initial_pop, temp_lhs_initial_pop)

        }

        if (verbose) {

            # updating the result_nrow variable
            result_nrow <- nrow(lhs_initial_pop)

            # printing progress
            print(result_nrow)

        }

    }

    # once we have the requested number of starting values, returning the final matrix
    return (lhs_initial_pop)

}
