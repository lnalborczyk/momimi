#' Models of motor inhibition during motor imagery
#'
#' Simulating data from the "threshold modulation model" (TMM) of motor imagery.
#'
#' @param nsims Numeric, number of simulations (observations/trials).
#' @param nsamples Numeric, number of samples (time steps) within a trial.
#' @param exec_threshold Numeric, motor execution threshold.
#' @param imag_threshold Numeric, motor imagery threshold.
#' @param peak_time Numeric, peak time of the activation function.
#' @param curvature Numeric, curvature of the activation function.
#' @param bw_noise Numeric, amount of between-trial noise.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param full_output Boolean, indicating whether activation curves should be returned.
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
#'     peak_time = log(0.5), curvature = 0.4, bw_noise = 0.08,
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
        nsamples = 5000,
        exec_threshold = 1,
        imag_threshold = 0.5,
        peak_time = 0,
        curvature = 0.6,
        bw_noise = NULL,
        uncertainty = c("par_level", "func_level", "diffusion"),
        full_output = FALSE,
        time_step = 0.001
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    # defining/retrieving the amount of between-trial noise (if any)
    bw_noise <- ifelse(
        test = uncertainty != "diffusion" & is.null(bw_noise),
        yes = 0.1, no = bw_noise
        )

    # if full_output = TRUE, returns the full activation function
    if (full_output == TRUE) {

        # computing the activation function and
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
                    peak_time = peak_time,
                    curvature = curvature,
                    uncertainty = uncertainty,
                    bw_noise = bw_noise,
                    time_step = time_step
                    )
                ) %>%
            # numerically finding the onset (RT) and offset
            dplyr::mutate(onset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::last() ) %>%
            dplyr::mutate(onset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::last() ) %>%
            # MT is defined as offset minus onset
            dplyr::mutate(mt_exec = .data$offset_exec - .data$onset_exec) %>%
            dplyr::mutate(mt_imag = .data$offset_imag - .data$onset_imag) %>%
            # convert from ms to seconds
            dplyr::mutate(dplyr::across(.data$onset_exec:.data$mt_imag, ~ . * time_step) ) %>%
            dplyr::ungroup()

        # setting the class of the resulting object
        class(results) <- c("momimi_full", "data.frame")

    } else if (full_output == FALSE) {

        if (uncertainty == "par_level") {

            # defining the activation rescaled lognormal function
            activation_function <- function (
                exec_threshold = 1, imag_threshold = 0.5,
                amplitude = 1.5, peak_time = 0,
                curvature = 0.4, bw_noise = 0.1
                ) {

                # adding some variability in the peak time and curvature parameters
                peak_time_sim <- stats::rnorm(n = 1, mean = peak_time, sd = bw_noise)
                curvature_sim <- stats::rnorm(n = 1, mean = curvature, sd = bw_noise)

                # no variability in the motor thresholds
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
                            peak_time = peak_time,
                            curvature = curvature,
                            bw_noise = bw_noise,
                            exec_threshold = .data$exec_threshold,
                            imag_threshold = .data$imag_threshold
                            )
                        )
                    ) %>%
                dplyr::ungroup()

            } else {

                # computing the activation function and
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
                            peak_time = peak_time,
                            curvature = curvature,
                            uncertainty = uncertainty,
                            bw_noise = bw_noise,
                            time_step = time_step
                            )
                        ) %>%
                    # numerically finding the onset (RT) and offset
                    dplyr::mutate(onset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::first() ) %>%
                    dplyr::mutate(offset_exec = which(.data$activation >= .data$exec_threshold) %>% dplyr::last() ) %>%
                    dplyr::mutate(onset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::first() ) %>%
                    dplyr::mutate(offset_imag = which(.data$activation >= .data$imag_threshold) %>% dplyr::last() ) %>%
                    # MT is defined as offset minus onset
                    dplyr::mutate(mt_exec = .data$offset_exec - .data$onset_exec) %>%
                    dplyr::mutate(mt_imag = .data$offset_imag - .data$onset_imag) %>%
                    # convert from ms to seconds
                    dplyr::mutate(dplyr::across(.data$onset_exec:.data$mt_imag, ~ . * time_step) ) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(.data$sim, .data$onset_imag, .data$mt_imag, .data$onset_exec, .data$mt_exec) %>%
                    dplyr::distinct()

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

    if (method == "functions") {

        x %>%
            tidyr::pivot_longer(., cols = .data$activation) %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = .data$time, y = .data$value,
                    group = interaction(.data$sim, .data$name)
                    )
                ) +
            # plotting the motor execution and motor imagery thresholds
            geomtextpath::geom_labelhline(
                yintercept = unique(x$exec_threshold),
                linetype = 2,
                hjust = 0.9,
                label = "Motor execution threshold"
                ) +
            geomtextpath::geom_labelhline(
                yintercept = unique(x$imag_threshold),
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
                colour = "black",
                linewidth = 1, alpha = 1,
                show.legend = FALSE
                ) +
            # ggplot2::ylim(c(0, 1.5) ) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                title = "Simulating activation patterns",
                x = "Time within a trial (s)",
                y = "Activation (a.u.)",
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
                x = "Reaction/Movement time (s)",
                y = "Probability density"
                )

    }

}

#' @export

plot.momimi <- function (x, ...) {

    plot.momimi_full(x = x, method = "distributions")

}
