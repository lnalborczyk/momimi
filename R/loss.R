#' Loss functions
#'
#' Provides possible loss functions for the two versions of the model.
#'
#' @param par Numeric, vector of parameter values.
#' @param data Dataframe, data to be used for computing the error.
#' @param nsims number of simulations (observations/trials).
#' @param nsamples number of samples (time steps) within a trial.
#' @param model_version Character, 3-par or 4-par threshold modulation model ("TMM3" or "TMM4").
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param time_step Numeric, time step used to numerical approximation.
#' @param diffusion_coef Numeric, diffusion coefficient.
#' @param exec_threshold motor execution threshold.
#' @param imag_threshold motor imagery threshold.
#' @param error_function Character, loss function to be used when fitting the model.
#'
#' @return The loss/error for a given set of observations.
#'
#' @importFrom rlang .data
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
        model_version = c("TMM3", "TMM4", "PIM"),
        exec_threshold = 1,
        imag_threshold = 0.5,
        uncertainty = c("par_level", "func_level", "diffusion"),
        time_step = 0.001,
        diffusion_coef = 0.001,
        error_function = c("g2", "rmse", "sse", "wsse", "ks")
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )

    # model_version should be one of above
    model_version <- match.arg(model_version)

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    # error_function should be one of above
    error_function <- match.arg(error_function)

    # how many trials should we simulate? if null, by default nrow(data)
    if (is.null(nsims) ) nsims <- as.numeric(nrow(data) )

    # what version of the model?
    if (model_version %in% c("TMM3", "TMM4") ) {

        if (model_version == "TMM3") {

            # retrieving parameter values for the activation function
            peak_time <- log(par[[2]])
            curvature <- par[[3]]

            # trying out the normalised lognormal distribution
            # amplitude_activ <- 1 / (peak_time_activ * curvature_activ * sqrt(2 * pi) )

            # retrieving parameter values for the execution threshold
            exec_threshold <- par[[1]]

            # defining imagery threshold relative to execution threshold
            imag_threshold <- imag_threshold * exec_threshold

            # defining the amount of between-trial noise
            bw_noise <- 0.1

        } else if (model_version == "TMM4") {

            # retrieving parameter values for the activation function
            peak_time <- log(par[[2]])
            curvature <- par[[3]]

            # retrieving parameter values for the execution threshold
            exec_threshold <- par[[1]]

            # defining imagery threshold relative to execution threshold
            imag_threshold <- imag_threshold * exec_threshold

            # retrieving the amount of between-trial variability
            bw_noise <- par[[4]]

        }

        ################################################################################
        # adding some constraints
        # ----------------------------------------------------------------------------
        # predicted RTs/MTs should be valid (not a NaN)
        # imagery threshold cannot be higher than execution threshold
        # balance max should not be above exec_threshold in imagined trials
        # balance max should not be above 4 * exec_threshold in executed trials
        # balance value at the end of the trial should be below 0.25 * exec_threshold
        ##############################################################################

        # defining the activation/inhibition rescaled lognormal function
        activation_function <- function (
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude = 1.5, peak_time = 0, curvature = 0.4,
        bw_noise = 0.1
        ) {

            # adding some variability in the parameters
            # amplitude_sim <- stats::rnorm(n = 1, mean = amplitude, sd = bw_noise)
            peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = bw_noise)
            curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = bw_noise)

            # no variability in the motor execution and motor imagery thresholds
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
        if (uncertainty == "par_level") {

            predicted_rt_mt <- suppressWarnings(activation_function(
                peak_time = peak_time,
                curvature = curvature,
                exec_threshold = exec_threshold,
                imag_threshold = imag_threshold,
                bw_noise = bw_noise
                ) )

        } else {

            # computing the activation/inhibition balance and
            # implied distributions of RTs and MTs per simulation
            predicted_rt_mt <- data.frame(
                time = 1:nsamples * time_step,
                exec_threshold = exec_threshold,
                imag_threshold = imag_threshold
                ) %>%
                dplyr::mutate(
                    activation = activation(
                        time = .data$time,
                        peak_time = peak_time,
                        curvature = curvature,
                        uncertainty = uncertainty,
                        bw_noise = bw_noise,
                        time_step = time_step,
                        diffusion_coef = diffusion_coef
                        )
                    ) %>%
                # numerically finding the balance's onset (RT) and offset
                dplyr::mutate(onset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::first() ) %>%
                dplyr::mutate(offset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::last() ) %>%
                # MT is defined as offset minus onset
                dplyr::mutate(mt_exec = .data$offset_exec - .data$onset_exec) %>%
                dplyr::mutate(onset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::first() ) %>%
                dplyr::mutate(offset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::last() ) %>%
                dplyr::mutate(mt_imag = .data$offset_imag - .data$onset_imag) %>%
                # convert from ms to seconds
                dplyr::mutate(dplyr::across(.data$onset_exec:.data$mt_imag, ~ . * time_step) ) %>%
                dplyr::select(.data$onset_imag, .data$mt_imag, .data$onset_exec, .data$mt_exec) %>%
                dplyr::distinct()

        }

        if (unique(data$action_mode) == "imagined") {

            predicted_rt_mt <- predicted_rt_mt[1:2]

        } else if (unique(data$action_mode) == "executed") {

            predicted_rt_mt <- predicted_rt_mt[3:4]

        }

        # coding the constraints (penalising by setting the prediction error to +Inf)
        if (any(is.na(predicted_rt_mt) ) ) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (imag_threshold >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        } else if (unique(data$action_mode) == "imagined" & !is.na(peak_time) & peak_time >= exec_threshold) {

            prediction_error <- Inf
            return (prediction_error)

        }

        # simulating some data
        # if uncertainty = "par_level", we can compute RT/MT analytically
        if (uncertainty == "par_level") {

            results <- data.frame(
                sim = rep(1:nsims, each = nsamples),
                exec_threshold = exec_threshold,
                imag_threshold = imag_threshold
                ) %>%
                dplyr::group_by(.data$sim) %>%
                dplyr::do(
                    suppressWarnings(
                        activation_function(
                            peak_time = peak_time,
                            curvature = curvature,
                            exec_threshold = .data$exec_threshold,
                            imag_threshold = .data$imag_threshold,
                            bw_noise = bw_noise
                            )
                        )
                    ) %>%
                dplyr::ungroup()

        # if uncertainty != "par_level", we need to compute RT/MT numerically
        } else {

            results <- data.frame(
                sim = rep(1:nsims, each = nsamples),
                sample = rep(1:nsamples, nsims),
                time = rep(1:nsamples, nsims) * time_step,
                exec_threshold = exec_threshold,
                imag_threshold = imag_threshold
                ) %>%
                dplyr::group_by(.data$sim) %>%
                dplyr::mutate(
                    activation = activation(
                        time = .data$time,
                        peak_time = peak_time,
                        curvature = curvature,
                        uncertainty = uncertainty,
                        bw_noise = bw_noise,
                        time_step = time_step,
                        diffusion_coef = diffusion_coef
                        )
                    ) %>%
                # numerically finding the balance's onset (RT) and offset
                dplyr::mutate(onset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::first() ) %>%
                dplyr::mutate(offset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::last() ) %>%
                # MT is defined as offset minus onset
                dplyr::mutate(mt_exec = .data$offset_exec - .data$onset_exec) %>%
                dplyr::mutate(onset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::first() ) %>%
                dplyr::mutate(offset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::last() ) %>%
                dplyr::mutate(mt_imag = .data$offset_imag - .data$onset_imag) %>%
                # convert from ms to seconds
                dplyr::mutate(dplyr::across(.data$onset_exec:.data$mt_imag, ~ . * time_step) ) %>%
                dplyr::ungroup() %>%
                dplyr::select(.data$sim, .data$onset_imag, .data$mt_imag, .data$onset_exec, .data$mt_exec) %>%
                dplyr::distinct()

        }

    }

    # retrieving distribution of simulated RTs
    if (unique(data$action_mode) == "imagined") {

        # retrieving distribution of simulated RTs
        predicted_rt <- tidyr::replace_na(data = results$onset_imag, replace = 1e6)

        # retrieving distribution of simulated MTs
        predicted_mt <- tidyr::replace_na(data = results$mt_imag, replace = 1e6)

    } else if (unique(data$action_mode) == "executed") {

        # retrieving distribution of simulated RTs
        predicted_rt <- tidyr::replace_na(data = results$onset_exec, replace = 1e6)

        # retrieving distribution of simulated MTs
        predicted_mt <- tidyr::replace_na(data = results$mt_exec, replace = 1e6)

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
        observed_rt_quantiles <- stats::quantile(x = data$reaction_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        observed_mt_quantiles <- stats::quantile(x = data$movement_time, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_rt_quantiles <- stats::quantile(x = predicted_rt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)
        predicted_mt_quantiles <- stats::quantile(x = predicted_mt, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), na.rm = TRUE)

        # quantile weights (cf. Ratcliff & Smith, 2004)
        quantile_weights <- c(2, 2, 1, 1, 0.5)

        # computing the weighted SSE
        prediction_error <- sum(quantile_weights * (predicted_rt_quantiles - observed_rt_quantiles)^2) +
            sum(quantile_weights * (predicted_mt_quantiles - observed_mt_quantiles)^2)

    } else if (error_function == "ks") {

        # computes the KS statistic
        prediction_error <- as.numeric(stats::ks.test(x = predicted_rt, y = data$reaction_time)$statistic) +
            as.numeric(stats::ks.test(x = predicted_mt, y = data$movement_time)$statistic)

    }

    # if NA, replacing it by +Inf
    if (is.na(prediction_error) ) prediction_error <- Inf

    # returning the prediction error
    return (prediction_error)

}
