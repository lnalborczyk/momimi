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
#'
#' @return ...
#'
#' @importFrom stats rnorm
#' @importFrom magrittr %>%
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

generating_initialpop <- function (
        nstudies, action_mode,
        par_names, lower_bounds, upper_bounds,
        model_version = c("TMM", "PIM")
        ) {

    # some tests for variable types
    stopifnot("nstudies must be a numeric..." = is.numeric(nstudies) )
    stopifnot("action_mode must be a character..." = is.character(action_mode) )

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

        # defining a function to compute the predicted RT and MT (quadratic formula)
        onset_offset <- function (alpha_f, alpha_g, mu_f, mu_g, sigma_f, sigma_g, thresh) {

            a <- sigma_g^2 - sigma_f^2
            b <- 2 * (sigma_f^2 * mu_g - sigma_g^2 * mu_f)
            c <- sigma_f^2 * mu_g^2 - sigma_g^2 * mu_f^2 - 2 * sigma_f^2 * sigma_g^2 * (log(alpha_f / alpha_g) - log(thresh) )
            onset <- exp((-b - sqrt(b^2 - 4 * a * c) ) / (2 * a) )
            offset <- exp((-b + sqrt(b^2 - 4 * a * c) ) / (2 * a) )

            return (c(onset, offset) )

        }

        # defining the balance function
        # basically a ratio of two rescaled lognormal functions
        balance_function <- function (
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
        amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6
        ) {

            # adding some variability in the other parameters
            # variability is currently fixed but could also be estimated
            amplitude_activ_sim <- rnorm(n = 1, mean = amplitude_activ, sd = 0.01)
            peak_time_activ_sim <- rnorm(n = 1, mean = peak_time_activ, sd = 0.01)
            curvature_activ_sim <- rnorm(n = 1, mean = curvature_activ, sd = 0.01)

            amplitude_inhib_sim <- rnorm(n = 1, mean = amplitude_inhib, sd = 0.01)
            peak_time_inhib_sim <- rnorm(n = 1, mean = peak_time_inhib, sd = 0.01)
            curvature_inhib_sim <- rnorm(n = 1, mean = curvature_inhib, sd = 0.01)

            # computing the predicted RT and MT in imagery
            onset_offset_imag <- onset_offset(
                alpha_f = amplitude_activ_sim, alpha_g = amplitude_inhib_sim,
                mu_f = peak_time_activ_sim, mu_g = peak_time_inhib_sim,
                sigma_f = curvature_activ_sim, sigma_g = curvature_inhib_sim,
                thresh = imag_threshold
                )

            onset_imag <- min(onset_offset_imag)
            mt_imag <- max(onset_offset_imag) - min(onset_offset_imag)

            # computing the predicted RT and MT in execution
            onset_offset_exec <- onset_offset(
                alpha_f = amplitude_activ_sim, alpha_g = amplitude_inhib_sim,
                mu_f = peak_time_activ_sim, mu_g = peak_time_inhib_sim,
                sigma_f = curvature_activ_sim, sigma_g = curvature_inhib_sim,
                thresh = exec_threshold
                )

            onset_exec <- min(onset_offset_exec)
            mt_exec <- max(onset_offset_exec) - min(onset_offset_exec)

            # returning it
            return (data.frame(onset_imag, mt_imag, onset_exec, mt_exec) )

        }

        # computing the predicted RT and MT
        predicted_rt_mt <- lhs_pars %>%
            rowwise() %>%
            do(
                suppressWarnings(
                    balance_function(
                        amplitude_activ = 1.5,
                        peak_time_activ = log(.$peak_time),
                        curvature_activ = .$curvature_activ,
                        amplitude_inhib = 1.5 / .$amplitude_ratio,
                        # peak_time_inhib = log(.$peak_time_inhib * .$peak_time_activ),
                        peak_time_inhib = log(.$peak_time),
                        curvature_inhib = .$curvature_inhib * .$curvature_activ
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
        # predicted RT should be be between 0.2 and 1 seconds
        # predicted MT should be be between 0.2 and 2 seconds
        # balance at the end of the trial should come back below 0.25
        final_par_values <- bind_cols(lhs_pars, predicted_rt_mt) %>%
            rowwise() %>%
            mutate(
                balance_end_of_trial = (amplitude_ratio) *
                    exp(-(log(3) - peak_time)^2 / (2 * curvature_activ^2) +
                            (log(3) - peak_time)^2 / (2 * curvature_inhib^2) )
            ) %>%
            # data.frame()
            mutate(
                included = case_when(
                    any(is.na(pick(everything() ) ) ) ~ FALSE,
                    pick(length(par_names) + 1) < 0.1 ~ FALSE,
                    pick(length(par_names) + 1) > 2 ~ FALSE,
                    pick(length(par_names) + 2) < 0.1 ~ FALSE,
                    pick(length(par_names) + 2) > 2 ~ FALSE,
                    balance_end_of_trial > 0.25 ~ FALSE,
                    .default = TRUE
                )
            )

        # reshaping the final matrix
        temp_lhs_initial_pop <- final_par_values %>%
            filter(included) %>%
            dplyr::select(1:length(par_names) ) %>%
            data.frame() %>%
            as.matrix()

        # rbinding results
        if (!exists("lhs_initial_pop") ) {

            lhs_initial_pop <- temp_lhs_initial_pop

        } else {

            lhs_initial_pop <- rbind(lhs_initial_pop, temp_lhs_initial_pop)

        }

        # updating the result_nrow variable
        result_nrow <- nrow(lhs_initial_pop)

        # printing progress
        print(result_nrow)

    }

    # once we have the requested number of starting values, returning the final matrix
    return (lhs_initial_pop)

}
