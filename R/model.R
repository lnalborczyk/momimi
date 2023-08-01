#' Models of motor inhibition during motor imagery
#'
#' Sampling from the two versions of the model.
#'
#' @param nsims Numeric, number of simulations (observations/trials).
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param exec_threshold Numeric, motor execution threshold.
#' @param imag_threshold Numeric, motor imagery threshold.
#' @param amplitude_activ Numeric, amplitude of the activation function.
#' @param peak_time_activ Numeric, peak time of the activation function.
#' @param curvature_activ Numeric, curvature of the activation function.
#' @param amplitude_inhib Numeric, amplitude of the inhibition function.
#' @param peak_time_inhib Numeric, peak time of the inhibition function.
#' @param curvature_inhib Numeric, curvature of the inhibition function.
#' @param bw_noise Numeric, amount of between-trial noise.
#' @param model_version Character, threshold modulation model ("TMM3" or "TMM4") or parallel inhibition model ("PIM").
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param full_output Boolean, indicating whether activ/inhib curves should be returned.
#' @param time_step Numeric, time step used to numerical approximation.
#'
#' @return A dataframe containing observation (i.e., RTs and MTs) and/or values of the underlying functions at each time step.
#'
#' @importFrom stats rnorm median
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # simulating 100 observations
#' simulated_data <- model(
#'     nsims = 100, nsamples = 2000,
#'     exec_threshold = 1, imag_threshold = 0.5,
#'     amplitude_activ = 0.8, peak_time_activ = log(0.5), curvature_activ = 0.4,
#'     amplitude_inhib = 1.75, peak_time_inhib = log(0.5), curvature_inhib = 0.5,
#'     model_version = "PIM",
#'     full_output = FALSE
#'     )
#'
#' # displaying the first six observations
#' head(simulated_data)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

model <- function (
        nsims = 100,
        nsamples = 3000,
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
        amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6,
        bw_noise = NULL,
        model_version = c("TMM3", "TMM4", "PIM"),
        uncertainty = c("par_level", "func_level", "diffusion"),
        full_output = FALSE,
        time_step = 0.001
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )
    stopifnot("curvature_inhib must be larger than curvature_activ..." = curvature_inhib > curvature_activ)

    # model_version should be one of above
    model_version <- match.arg(model_version)

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    # defining/retrieving the amount of between-trial noise (if any)
    bw_noise <- ifelse(
        test = uncertainty != "diffusion" & is.null(bw_noise),
        yes = 0.1, no = bw_noise
        )

    # if full_output = TRUE, returns the full activation, inhibition,
    # and balance functions
    if (full_output == TRUE) {

        # computing the activation/inhibition balance and
        # implied distributions of RTs and MTs per simulation
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
                    amplitude = amplitude_activ,
                    peak_time = peak_time_activ,
                    curvature = curvature_activ,
                    uncertainty = uncertainty,
                    bw_noise = bw_noise,
                    time_step = time_step
                    )
                ) %>%
            dplyr::mutate(
                inhibition = activation(
                    time = .data$time,
                    amplitude = amplitude_inhib,
                    peak_time = peak_time_inhib,
                    curvature = curvature_inhib,
                    uncertainty = uncertainty,
                    bw_noise = bw_noise,
                    time_step = time_step
                    )
                ) %>%
            {if (model_version == "PIM") dplyr::mutate(., balance = .data$activation / .data$inhibition) else dplyr::mutate(., balance = .data$activation)} %>%
            # numerically finding the balance's onset (RT) and offset
            dplyr::mutate(onset_exec = which(.data$balance >= .data$exec_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_exec = which(.data$balance >= .data$exec_threshold) %>% dplyr::last() ) %>%
            # MT is defined as offset minus onset
            dplyr::mutate(mt_exec = .data$offset_exec - .data$onset_exec) %>%
            dplyr::mutate(onset_imag = which(.data$balance >= .data$imag_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_imag = which(.data$balance >= .data$imag_threshold) %>% dplyr::last() ) %>%
            dplyr::mutate(mt_imag = .data$offset_imag - .data$onset_imag) %>%
            # convert from ms to seconds
            dplyr::mutate(dplyr::across(.data$onset_exec:.data$mt_imag, ~ . * time_step) ) %>%
            dplyr::ungroup()

        # setting the class of the resulting object
        class(results) <- c("momimi_full", "data.frame")

    } else if (full_output == FALSE) {

        if (model_version %in% c("TMM3", "TMM4") ) {

            if (uncertainty == "par_level") {

                # defining the activation/inhibition rescaled lognormal function
                activation_function <- function (
                    exec_threshold = 1, imag_threshold = 0.5,
                    amplitude = 1.5, peak_time = 0,
                    curvature = 0.4, bw_noise = 0.1
                    ) {

                    # adding some variability in the other parameters
                    # variability is currently fixed but could also be estimated
                    amplitude_sim <- stats::rnorm(n = 1, mean = amplitude, sd = bw_noise)
                    peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = bw_noise)
                    curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = bw_noise)

                    # no variability in the motor imagery threshold
                    exec_threshold_sim <- exec_threshold
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

                # computing the activation function
                # and the implied/predicted distributions of RTs and MTs
                results <- data.frame(
                    sim = rep(1:nsims, each = nsamples),
                    exec_threshold = exec_threshold,
                    imag_threshold = imag_threshold
                    ) %>%
                    dplyr::group_by(.data$sim) %>%
                    dplyr::do(
                        suppressWarnings(
                            activation_function(
                                amplitude = amplitude_activ,
                                peak_time = peak_time_activ,
                                curvature = curvature_activ,
                                bw_noise = bw_noise,
                                exec_threshold = .data$exec_threshold,
                                imag_threshold = .data$imag_threshold
                                )
                            )
                        ) %>%
                    dplyr::ungroup()

                } else {

                    # computing the activation/inhibition balance and
                    # implied distributions of RTs and MTs per simulation
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
                                amplitude = amplitude_activ,
                                peak_time = peak_time_activ,
                                curvature = curvature_activ,
                                uncertainty = uncertainty,
                                bw_noise = bw_noise,
                                time_step = time_step
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

        } else if (model_version == "PIM") {

            # defining the balance function
            # basically a ratio of two rescaled lognormal functions
            balance_function <- function (
            exec_threshold = 1, imag_threshold = 0.5,
            amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
            amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6,
            bw_noise = 0.1
            ) {

                # adding some variability in the other parameters
                # variability is currently fixed but could also be estimated
                amplitude_activ_sim <- stats::rnorm(n = 1, mean = amplitude_activ, sd = bw_noise)
                peak_time_activ_sim <- stats::rnorm(n = 1, mean = peak_time_activ, sd = bw_noise)
                curvature_activ_sim <- stats::rnorm(n = 1, mean = curvature_activ, sd = bw_noise)

                amplitude_inhib_sim <- stats::rnorm(n = 1, mean = amplitude_inhib, sd = bw_noise)
                peak_time_inhib_sim <- stats::rnorm(n = 1, mean = peak_time_inhib, sd = bw_noise)
                curvature_inhib_sim <- stats::rnorm(n = 1, mean = curvature_inhib, sd = bw_noise)

                # computing the predicted RT and MT in imagery
                onset_offset_imag <- onset_offset(
                    amplitude_activ = amplitude_activ_sim,
                    amplitude_inhib = amplitude_inhib_sim,
                    peak_time_activ = peak_time_activ_sim,
                    peak_time_inhib = peak_time_inhib_sim,
                    curvature_activ = curvature_activ_sim,
                    curvature_inhib = curvature_inhib_sim,
                    thresh = imag_threshold,
                    model_version = model_version
                    )

                onset_imag <- min(onset_offset_imag)
                mt_imag <- max(onset_offset_imag) - min(onset_offset_imag)

                # computing the predicted RT and MT in execution
                onset_offset_exec <- onset_offset(
                    amplitude_activ = amplitude_activ_sim,
                    amplitude_inhib = amplitude_inhib_sim,
                    peak_time_activ = peak_time_activ_sim,
                    peak_time_inhib = peak_time_inhib_sim,
                    curvature_activ = curvature_activ_sim,
                    curvature_inhib = curvature_inhib_sim,
                    thresh = exec_threshold,
                    model_version = model_version
                    )

                onset_exec <- min(onset_offset_exec)
                mt_exec <- max(onset_offset_exec) - min(onset_offset_exec)

                # returning it
                return (data.frame(onset_imag, mt_imag, onset_exec, mt_exec) )

            }

            # computing the activation/inhibition balance
            # and the implied/predicted distributions of RTs and MTs
            results <- data.frame(
                sim = rep(1:nsims, each = nsamples),
                exec_threshold = exec_threshold,
                imag_threshold = imag_threshold
                ) %>%
                dplyr::group_by(.data$sim) %>%
                dplyr::do(
                    suppressWarnings(
                        balance_function(
                            amplitude_activ = amplitude_activ,
                            peak_time_activ = peak_time_activ,
                            curvature_activ = curvature_activ,
                            amplitude_inhib = amplitude_inhib,
                            peak_time_inhib = peak_time_inhib,
                            curvature_inhib = curvature_inhib,
                            exec_threshold = .data$exec_threshold,
                            imag_threshold = .data$imag_threshold
                            )
                        )
                    ) %>%
                dplyr::ungroup()

        }

        # setting the class of the resulting object
        class(results) <- c("momimi", "data.frame")

    }

    # returning the results
    return (results)

}

#' @export

plot.momimi_full <- function (x, method = c("functions", "distributions"), ...) {

    # ensuring that method is one of the above
    method <- match.arg(method)

    # testing whether activation is different than balance (i.e., TMM or PIM?)
    multiple_functions <- any(x$activation != x$balance)

    if (method == "functions") {

        x %>%
            {if (any(x$activation == x$balance) ) tidyr::pivot_longer(., cols = .data$activation) else tidyr::pivot_longer(., cols = .data$activation:.data$balance)} %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = .data$time, y = .data$value,
                    group = interaction(.data$sim, .data$name),
                    colour = .data$name
                    )
                ) +
            # plotting the motor execution and motor imagery thresholds
            geomtextpath::geom_labelhline(
                yintercept = 1,
                linetype = 2,
                hjust = 0.9,
                label = "Motor execution threshold"
                ) +
            geomtextpath::geom_labelhline(
                yintercept = 0.5,
                linetype = 2,
                hjust = 0.9,
                label = "Motor imagery threshold"
                ) +
            # plotting some individual simulations
            ggplot2::geom_line(
                data = . %>% dplyr::filter(.data$sim %in% unique(.data$sim)[1:20]),
                linewidth = 0.5, alpha = 0.2,
                # colour = "grey",
                show.legend = FALSE
                ) +
            # plotting average
            ggplot2::stat_summary(
                ggplot2::aes(group = .data$name, colour = .data$name),
                fun = "median", geom = "line",
                # colour = "black",
                linewidth = 1, alpha = 1,
                show.legend = multiple_functions
                ) +
            # ggplot2::ylim(c(0, 1.5) ) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                title = "Simulating activation/inhibition patterns",
                x = "Time within a trial (in seconds)",
                y = "Activation/inhibition (a.u.)",
                colour = "",
                fill = ""
                )

    } else if (method == "distributions") {

        x %>%
            dplyr::mutate(
                exec_rt_median = stats::median(.data$onset_exec),
                imag_rt_median = stats::median(.data$onset_imag),
                exec_mt_median = stats::median(.data$mt_exec),
                imag_mt_median = stats::median(.data$mt_imag)
                ) %>%
            tidyr::pivot_longer(cols = c(.data$onset_imag, .data$mt_imag) ) %>%
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
            ggplot2::geom_label(
                data = . %>% dplyr::summarise(m = unique(.data$imag_rt_median) ),
                ggplot2::aes(x = .data$m, y = 0, label = round(.data$m, 3) ),
                position = ggplot2::position_nudge(y = 0.01),
                size = 4,
                inherit.aes = FALSE
                ) +
            ggplot2::geom_label(
                data = . %>% dplyr::summarise(m = unique(.data$imag_mt_median) ),
                ggplot2::aes(x = .data$m, y = 0, label = round(.data$m, 3) ),
                position = ggplot2::position_nudge(y = 0.01),
                size = 4,
                inherit.aes = FALSE
                ) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                title = "Simulating the implied distributions of RTs and MTs",
                x = "Reaction/Movement time (in seconds)",
                y = "Probability density"
                )

    }

}

#' @export

plot.momimi <- function (x, ...) {

    plot.momimi_full(x = x, method = "distributions")

}
