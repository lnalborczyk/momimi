#' Generating initial population
#'
#' Generating plausible initial parameter values for DEoptim based on empirical constraints.
#'
#' @param nstudies Numeric, number of starting values in the LHS.
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param action_mode Character, action mode (executed or imagined).
#' @param par_names Character, vector of parameter names.
#' @param lower_bounds Numeric, vector of lower bounds for parameters.
#' @param upper_bounds Numeric, vector of upper bounds for parameters.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param time_step Numeric, time step used to numerical approximation.
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
        nstudies,
        nsamples = 5000,
        action_mode,
        par_names, lower_bounds, upper_bounds,
        uncertainty = c("par_level", "func_level", "diffusion"),
        time_step = 0.001,
        rt_contraints = c(0.1, 2),
        mt_contraints = c(0.1, 2),
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

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

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

        # defining the activation/inhibition rescaled lognormal function
        activation_function <- function (
            exec_threshold = 1, imag_threshold = 0.5,
            amplitude = 1.5, peak_time = 0, curvature = 0.4, bw_noise = 0.1
            ) {

            # adding some variability in the other parameters
            # variability is currently fixed but could also be estimated
            peak_time_sim <- rnorm(n = 1, mean = peak_time, sd = bw_noise)
            curvature_sim <- rnorm(n = 1, mean = curvature, sd = bw_noise)

            # # no variability in the motor execution or motor imagery thresholds
            exec_threshold_sim <- exec_threshold
            imag_threshold_sim <- imag_threshold

            # computing the predicted RT and MT in imagery
            onset_offset_imag <- onset_offset(
                peak_time = peak_time_sim,
                curvature = curvature_sim,
                thresh = imag_threshold_sim
                )

            onset_imag <- min(onset_offset_imag)
            mt_imag <- max(onset_offset_imag) - min(onset_offset_imag)

            # computing the predicted RT and MT in execution
            # computing the predicted RT and MT in imagery
            onset_offset_exec <- onset_offset(
                peak_time = peak_time_sim,
                curvature = curvature_sim,
                thresh = exec_threshold_sim
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
                        exec_threshold = .data$exec_threshold,
                        imag_threshold = 0.5 * .data$exec_threshold,
                        peak_time = log(.data$peak_time),
                        curvature = .data$curvature,
                        bw_noise = .data$bw_noise
                        )
                    )
                )

        if (action_mode == "imagined") {

            predicted_rt_mt <- predicted_rt_mt[, 1:2]

        } else if (action_mode == "executed") {

            predicted_rt_mt <- predicted_rt_mt[, 3:4]

        }

        # adding some constraints
        # to restrict parameter values to realistic RT/MT data (e.g., Evans, 2020)
        # predicted RT/MT should be valid (not a NaN)
        # predicted RT should be be between min(rt_constraints) and max(rt_constraints
        # predicted MT should be be between min(mt_constraints) and max(mt_constraints)
        final_par_values <- dplyr::bind_cols(lhs_pars, predicted_rt_mt) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                included = dplyr::case_when(
                    any(is.na(dplyr::pick(dplyr::everything() ) ) ) ~ FALSE,
                    dplyr::pick(length(par_names) + 1) < min(rt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 1) > max(rt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 2) < min(mt_contraints) ~ FALSE,
                    dplyr::pick(length(par_names) + 2) > max(mt_contraints) ~ FALSE,
                    .default = TRUE
                    )
                )

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
