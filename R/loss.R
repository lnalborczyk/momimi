#' Loss functions
#'
#' Provides possible loss functions for the two versions of the model.
#'
#' @param par Numeric, vector of parameter values.
#' @param data Dataframe, data to be used for computing the error.
#' @param nsims number of simulations (observations/trials).
#' @param nsamples number of samples (time steps) within a trial.
#' @param model_version Version of the model ("TMM" or "PIM").
#' @param exec_threshold motor execution threshold.
#' @param imag_threshold motor imagery threshold.
#' @param error_function Character, loss function to be used when fitting the model.
#'
#' @return A dataframe
#'
#' @importFrom magrittr %>%
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

# simulating some data and computing the prediction error
loss <- function (
        par, data,
        nsims = NULL,
        nsamples = 3000,
        model_version = "TMM",
        exec_threshold = 1, imag_threshold = 0.5,
        error_function = c("g2", "rmse", "sse", "wsse", "ks")
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )

    # error_function should be one of above
    error_function <- match.arg(error_function)

    # how many trials should we simulate? if null, by default nrow(data)
    if (is.null(nsims) ) nsims <- as.numeric(nrow(data) )

    # what version of the model?
    if (model_version == "TMM") {

        # setting an arbitrary value for the amplitude of the activation function
        amplitude_activ <- par[[1]]

        # retrieving parameter values for the activation function
        peak_time_activ <- log(par[[2]])
        curvature_activ <- par[[3]]

        # retrieving parameter values for the execution threshold
        exec_threshold <- par[[4]] * amplitude_activ

        # defining imagery threshold relative to execution threshold
        imag_threshold <- imag_threshold * exec_threshold

        ################################################################################
        # adding some constraints
        # ----------------------------------------------------------------------------
        # predicted RTs/MTs should be valid (not a NaN)
        # imagery threshold cannot be higher than execution threshold
        # balance max should not be above exec_threshold in imagined trials
        # balance max should not be above 4 * exec_threshold in executed trials
        # balance value at the end of the trial should be below 0.25
        ##########################################################################

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

            # no variability in the motor imagery threshold
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
        predicted_rt_mt <- suppressWarnings(activation_function(
            amplitude = amplitude_activ,
            peak_time = peak_time_activ,
            curvature = curvature_activ,
            exec_threshold = exec_threshold,
            imag_threshold = imag_threshold
            ) )

        if (unique(data$action_mode) == "imagined") {

            predicted_rt_mt <- predicted_rt_mt[1:2]

        } else if (unique(data$action_mode) == "executed") {

            predicted_rt_mt <- predicted_rt_mt[3:4]

        }

        # computing the balance output at the end of the trial (i.e., when time = 3)
        balance_end_of_trial <- amplitude_activ * exp(-(log(nsamples / 1e3) - peak_time_activ)^2 / (2 * curvature_activ^2) )

        # coding the constraints (penalising by setting the prediction error to +Inf)
        if (any(is.na(predicted_rt_mt) ) ) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (imag_threshold >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (unique(data$action_mode) == "imagined" & !is.na(peak_time_activ) & peak_time_activ >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (unique(data$action_mode) == "executed" & !is.na(peak_time_activ) & peak_time_activ >= 4 * exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (balance_end_of_trial >= 0.5 * imag_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        }

        # simulating some data
        results <- data.frame(
            sim = rep(1:nsims, each = nsamples),
            exec_threshold = exec_threshold,
            imag_threshold = imag_threshold
            ) %>%
            dplyr::group_by(sim) %>%
            dplyr::do(
                suppressWarnings(
                    activation_function(
                        amplitude = amplitude_activ,
                        peak_time = peak_time_activ,
                        curvature = curvature_activ,
                        exec_threshold = exec_threshold,
                        imag_threshold = imag_threshold
                        )
                    )
                ) %>%
            dplyr::ungroup()

    } else if (model_version == "PIM") {

        # defines imagery threshold relative to execution threshold
        imag_threshold <- imag_threshold * exec_threshold

        # setting an arbitrary value for the amplitude of the activation function
        amplitude_activ <- 1.5

        # retrieving parameter values for the activation function
        peak_time_activ <- log(par[[2]])
        curvature_activ <- par[[3]]

        # setting a value for the ratio amplitude_activ / amplitude_inhib
        amplitude_inhib <- amplitude_activ / par[[1]]

        # retrieving parameter values for the inhibition function
        peak_time_inhib <- log(par[[2]])
        curvature_inhib <- par[[4]] * par[[3]]

        ############################################################################
        # adding some constraints
        # -------------------------------------------------------------------------
        # predicted RTs/MTs should be valid (not a NaN)
        # curvature_activ should be lower than curvature_inhib
        # imagery threshold cannot be higher than execution threshold
        # balance max should not be above exec_threshold in imagined trials
        # balance max should not be above 4 * exec_threshold in executed trials
        # balance value at the end of the trial should be below 0.25
        ##########################################################################

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

            # in this model, there is no variation in the thresholds
            # exec_threshold_sim <- exec_threshold
            # imag_threshold_sim <- imag_threshold

            # computing the predicted RT and MT in imagery
            onset_offset_imag <- onset_offset(
                amplitude_activ = amplitude_activ_sim,
                peak_time_activ = peak_time_activ_sim,
                curvature_activ = curvature_activ_sim,
                amplitude_inhib = curvature_inhib_sim,
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
                curvature_activ = curvature_activ_sim,
                amplitude_inhib = curvature_inhib_sim,
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
        predicted_rt_mt <- suppressWarnings(balance_function(
            exec_threshold = exec_threshold,
            imag_threshold = imag_threshold,
            amplitude_activ = amplitude_activ,
            peak_time_activ = peak_time_activ,
            curvature_activ = curvature_activ,
            amplitude_inhib = amplitude_inhib,
            peak_time_inhib = peak_time_inhib,
            curvature_inhib = curvature_inhib
            ) )

        if (unique(data$action_mode) == "imagined") {

            predicted_rt_mt <- predicted_rt_mt[1:2]

        } else if (unique(data$action_mode) == "executed") {

            predicted_rt_mt <- predicted_rt_mt[3:4]

        }

        # computing the peak time (mode) of the balance function
        balance_peak_time <- exp(
            (peak_time_activ * curvature_inhib^2 - peak_time_inhib * curvature_activ^2) /
                (curvature_inhib^2 - curvature_activ^2) )

        # computing the maximum value of the balance function
        balance_max <- (amplitude_activ / amplitude_inhib) *
            exp(-(log(balance_peak_time) - peak_time_activ)^2 / (2 * curvature_activ^2) +
                    (log(balance_peak_time) - peak_time_inhib)^2 / (2 * curvature_inhib^2) )

        # computing the balance output at the end of the trial (i.e., when time = 3)
        balance_end_of_trial <- (amplitude_activ / amplitude_inhib) *
            exp(-(log(nsamples / 1e3) - peak_time_activ)^2 / (2 * curvature_activ^2) +
                    (log(nsamples / 1e3) - peak_time_inhib)^2 / (2 * curvature_inhib^2) )

        # coding the constraints (penalising by setting the prediction error to +Inf)
        if (any(is.na(predicted_rt_mt) ) ) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (curvature_activ >= curvature_inhib) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (imag_threshold >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (unique(data$action_mode) == "imagined" & !is.na(balance_max) & balance_max >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (unique(data$action_mode) == "executed" & !is.na(balance_max) & balance_max >= 4 * exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (balance_end_of_trial >= 0.25) {

            prediction_error <- Inf
            return (prediction_error)

        }

        # simulating some data
        results <- data.frame(
            sim = rep(1:nsims, each = nsamples),
            exec_threshold = exec_threshold,
            imag_threshold = imag_threshold
            ) |>
            dplyr::group_by(sim) |>
            dplyr::do(
                suppressWarnings(
                    balance_function(
                        exec_threshold = exec_threshold,
                        imag_threshold = imag_threshold,
                        amplitude_activ = amplitude_activ,
                        peak_time_activ = peak_time_activ,
                        curvature_activ = curvature_activ,
                        amplitude_inhib = amplitude_inhib,
                        peak_time_inhib = peak_time_inhib,
                        curvature_inhib = curvature_inhib
                        )
                    )
                ) |>
            dplyr::ungroup()

    }

    # retrieving distribution of simulated RTs
    if (unique(data$action_mode) == "imagined") {

        # retrieving distribution of simulated RTs
        predicted_rt <- replace_na(data = results$onset_imag, replace = 1e6)

        # retrieving distribution of simulated MTs
        predicted_mt <- replace_na(data = results$mt_imag, replace = 1e6)

    } else if (unique(data$action_mode) == "executed") {

        # retrieving distribution of simulated RTs
        predicted_rt <- replace_na(data = results$onset_exec, replace = 1e6)

        # retrieving distribution of simulated MTs
        predicted_mt <- replace_na(data = results$mt_exec, replace = 1e6)

    } else {

        warning ("Action-mode should be one of 'imagined' or 'executed'...")

    }

    if (error_function == "g2") {

        # what quantiles should we look at?
        quantile_probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)

        # computes observed RT quantiles
        observed_rt_quantiles <- stats::quantile(x = data$reaction_time, probs = quantile_probs, na.rm = TRUE)

        # computes observed MT quantiles
        observed_mt_quantiles <- stats::quantile(x = data$movement_time, probs = quantile_probs, na.rm = TRUE)

        # computes observed proportion of data in RT quantiles
        observed_rt_quantiles_props <- quantiles_props(
            x = data$reaction_time, quants = observed_rt_quantiles
            )

        # computes observed proportion of data in MT quantiles
        observed_mt_quantiles_props <- quantiles_props(
            x = data$movement_time, quants = observed_mt_quantiles
            )

        # computes predicted proportion of data in RT quantiles
        predicted_rt_quantiles_props <- quantiles_props(
            x = predicted_rt, quants = observed_rt_quantiles
            )

        # computes predicted proportion of data in MT quantiles
        predicted_mt_quantiles_props <- quantiles_props(
            x = predicted_mt, quants = observed_mt_quantiles
            )

        # applies a small correction when prop = 0 to avoid negative or Inf g-square
        predicted_rt_quantiles_props <- ifelse(
            test = predicted_rt_quantiles_props == 0,
            yes = 1e-6,
            no = predicted_rt_quantiles_props
            )

        # makes sure proportions sum to 1
        predicted_rt_quantiles_props <- predicted_rt_quantiles_props /
            sum(predicted_rt_quantiles_props)

        # applies a small correction when prop = 0 to avoid negative or Inf g-square
        predicted_mt_quantiles_props <- ifelse(
            test = predicted_mt_quantiles_props == 0,
            yes = 1e-6,
            no = predicted_mt_quantiles_props
            )

        # makes sure proportions sum to 1
        predicted_mt_quantiles_props <- predicted_mt_quantiles_props /
            sum(predicted_mt_quantiles_props)

        # computes the G^2 prediction error
        # which is the error for RTs plus the error for MTs
        # see Ratcliff & Smith (2004, doi:10.1037/0033-295X.111.2.333) or
        # Servant et al. (2019, doi:10.1152/jn.00507.2018)
        prediction_error <- 2 * sum(observed_rt_quantiles_props * log(observed_rt_quantiles_props / predicted_rt_quantiles_props) ) +
            sum(observed_mt_quantiles_props * log(observed_mt_quantiles_props / predicted_mt_quantiles_props) )

    } else if (error_function == "rmse") {

        # or RMSE as in Ulrich et al. (2016)
        observed_rt_quantiles <- stats::quantile(x = data$reaction_time, probs = seq(0, 1, 0.1), na.rm = TRUE)
        observed_mt_quantiles <- stats::quantile(x = data$movement_time, probs = seq(0, 1, 0.1), na.rm = TRUE)
        predicted_rt_quantiles <- stats::quantile(x = predicted_rt, probs = seq(0, 1, 0.1), na.rm = TRUE)
        predicted_mt_quantiles <- stats::quantile(x = predicted_mt, probs = seq(0, 1, 0.1), na.rm = TRUE)

        # computing the weighted RMSE
        prediction_error <- sqrt((sum(predicted_rt_quantiles - observed_rt_quantiles)^2 +
                                      sum(predicted_mt_quantiles - observed_mt_quantiles)^2) /
                                     (length(predicted_rt_quantiles) + length(predicted_mt_quantiles) ) )

    } else if (error_function == "sse") {

        # or raw SSE as in Ractliff & Smith (2004)
        observed_rt_quantiles <- stats::quantile(x = data$reaction_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        observed_mt_quantiles <- stats::quantile(x = data$movement_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_rt_quantiles <- stats::quantile(x = predicted_rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_mt_quantiles <- stats::quantile(x = predicted_mt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)

        # computing the weighted SSE
        prediction_error <- sum((predicted_rt_quantiles - observed_rt_quantiles)^2) +
            sum((predicted_mt_quantiles - observed_mt_quantiles)^2)

    } else if (error_function == "wsse") {

        # or weighted SSE as in Ratcliff & Smith (2004)
        observed_rt_quantiles <- quantile(x = data$reaction_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        observed_mt_quantiles <- quantile(x = data$movement_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_rt_quantiles <- quantile(x = predicted_rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_mt_quantiles <- quantile(x = predicted_mt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)

        # quantile weights (cf. Ratcliff & Smith, 2004)
        quantile_weights <- c(2, 2, 1, 1, 0.5)

        # computing the weighted SSE
        prediction_error <- sum(quantile_weights * (predicted_rt_quantiles - observed_rt_quantiles)^2) +
            sum(quantile_weights * (predicted_mt_quantiles - observed_mt_quantiles)^2)

    } else if (error_function == "ks") {

        # computes the KS statistic
        prediction_error <- as.numeric(ks.test(x = predicted_rt, y = data$reaction_time)$statistic) +
            as.numeric(ks.test(x = predicted_mt, y = data$movement_time)$statistic)

    }

    # if NA, replacing it by +Inf
    if (is.na(prediction_error) ) prediction_error <- Inf

    # returning the prediction error
    return (prediction_error)

}