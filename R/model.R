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
#' @param model_version Character, threshold modulation model ("TMM") or parallel inhibition model ("PIM").
#' @param full_output Boolean, indicating whether activ/inhib curves should be returned.
#'
#' @return A dataframe
#'
#' @importFrom stats rnorm median
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' simulation_results <- model(
#'     nsims = 100, nsamples = 2000,
#'     exec_threshold = 1, imag_threshold = 0.5,
#'     amplitude_activ = 0.8, peak_time_activ = log(0.5), curvature_activ = 0.4,
#'     amplitude_inhib = 1.75, peak_time_inhib = log(0.5), curvature_inhib = 0.5,
#'     model_version = "PIM",
#'     full_output = FALSE
#'     )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

model <- function (
        nsims = 100, nsamples = 3000,
        exec_threshold = 1, imag_threshold = 0.5,
        amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
        amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6,
        model_version = c("TMM", "PIM"),
        full_output = FALSE
        ) {

    # some tests for variable types
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nsamples must be a numeric..." = is.numeric(nsamples) )

    # model_version should be one of above
    model_version <- match.arg(model_version)

    # if full_output = TRUE, returns the full activation, inhibition,
    # and balance functions
    if (full_output == TRUE) {

        # computing the activation/inhibition balance and
        # implied distributions of RTs and MTs per simulation
        results <- data.frame(
            sim = rep(1:nsims, each = nsamples),
            sample = rep(1:nsamples, nsims),
            time = rep(1:nsamples, nsims) / 1e3,
            exec_threshold = exec_threshold,
            imag_threshold = imag_threshold
            ) %>%
            dplyr::group_by(sim) %>%
            dplyr::mutate(
                activation = activation(
                    time = time,
                    amplitude = amplitude_activ,
                    peak_time = peak_time_activ,
                    curvature = curvature_activ
                    )
                ) %>%
            dplyr::mutate(
                inhibition = activation(
                    time = time,
                    amplitude = amplitude_inhib,
                    peak_time = peak_time_inhib,
                    curvature = curvature_inhib
                    )
                ) %>%
            {if (model_version == "PIM") dplyr::mutate(., balance = activation / inhibition) else dplyr::mutate(., balance = activation)} %>%
            # numerically finding the balance's onset (RT) and offset
            dplyr::mutate(onset_exec = which(balance > exec_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_exec = which(balance > exec_threshold) %>% dplyr::last() ) %>%
            # MT is defined as offset minus onset
            dplyr::mutate(mt_exec = offset_exec - onset_exec) %>%
            dplyr::mutate(onset_imag = which(balance > imag_threshold) %>% dplyr::first() ) %>%
            dplyr::mutate(offset_imag = which(balance > imag_threshold) %>% dplyr::last() ) %>%
            dplyr::mutate(mt_imag = offset_imag - onset_imag) %>%
            # convert from ms to seconds
            dplyr::mutate(across(onset_exec:mt_imag, ~ . / 1e3) ) %>%
            dplyr::ungroup()

        # setting the class of the resulting object
        class(results) <- c("momimi_full", "data.frame")

    } else if (full_output == FALSE) {

        if (model_version == "TMM") {

            # defining a function to compute the predicted RT and MT (quadratic formula)
            # onset_offset <- function (alpha, mu, sigma, thresh) {
            #
            #     onset <-  exp(mu - sqrt(-2 * sigma^2 * log(thresh / alpha) ) )
            #     offset <- exp(mu + sqrt(-2 * sigma^2 * log(thresh / alpha) ) )
            #
            #     return (c(onset, offset) )
            #
            # }

            # defining the activation/inhibition rescaled lognormal function
            activation_function <- function (exec_threshold = 1,
                                             imag_threshold = 0.5,
                                             amplitude = 1.5, peak_time = 0,
                                             curvature = 0.4
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
                ungroup()

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
                amplitude_activ_sim <- rnorm(n = 1, mean = amplitude_activ, sd = 0.01)
                peak_time_activ_sim <- rnorm(n = 1, mean = peak_time_activ, sd = 0.01)
                curvature_activ_sim <- rnorm(n = 1, mean = curvature_activ, sd = 0.01)

                amplitude_inhib_sim <- rnorm(n = 1, mean = amplitude_inhib, sd = 0.01)
                peak_time_inhib_sim <- rnorm(n = 1, mean = peak_time_inhib, sd = 0.01)
                curvature_inhib_sim <- rnorm(n = 1, mean = curvature_inhib, sd = 0.01)

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
                dplyr::group_by(sim) %>%
                dplyr::do(
                    suppressWarnings(
                        balance_function(
                            amplitude_activ = amplitude_activ,
                            peak_time_activ = peak_time_activ,
                            curvature_activ = curvature_activ,
                            amplitude_inhib = amplitude_inhib,
                            peak_time_inhib = peak_time_inhib,
                            curvature_inhib = curvature_inhib,
                            exec_threshold = exec_threshold,
                            imag_threshold = imag_threshold
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

plot.momimi_full <- function (x, method = c("functions", "distributions", "both"), ...) {

    # ensuring that the method is one of the above
    method <- match.arg(method)

    # plotting the functions
    p1 <- x %>%
        # pivot_longer(cols = activation:balance) %>%
        tidyr::pivot_longer(cols = activation) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = time, y = value,
                group = interaction(sim, name),
                colour = name
                )
            ) +
        # plotting the motor execution and motor imagery thresholds
        geomtextpath::geom_labelhline(
            yintercept = 1, linetype = 2,
            hjust = 0.9,
            label = "Motor execution threshold"
            ) +
        geomtextpath::geom_labelhline(
            yintercept = 0.5, linetype = 2,
            hjust = 0.9,
            label = "Motor imagery threshold"
            ) +
        # plotting some individual simulations
        ggplot2::geom_line(
            data = . %>% dplyr::filter(sim %in% unique(sim)[1:50]),
            linewidth = 0.5, alpha = 0.5, colour = "grey",
            show.legend = FALSE
            ) +
        # plotting average
        ggplot2::stat_summary(
            ggplot2::aes(group = name, colour = name),
            fun = "median", geom = "line",
            colour = "black",
            linewidth = 1, alpha = 1,
            show.legend = TRUE
            ) +
        ggplot2::ylim(c(0, 1.5) ) +
        ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
        ggplot2::labs(
            # title = "Simulating activation/inhibition patterns",
            title = "Simulating activation patterns",
            x = "Time within a trial (in seconds)",
            # y = "Activation/inhibition (a.u.)",
            y = "Activation (a.u.)",
            colour = "",
            fill = ""
            )

    p2 <- x %>%
        dplyr::mutate(
            exec_rt_median = stats::median(onset_exec),
            imag_rt_median = stats::median(onset_imag),
            exec_mt_median = stats::median(mt_exec),
            imag_mt_median = stats::median(mt_imag),
            ) %>%
        tidyr::pivot_longer(cols = c(onset_imag, mt_imag) ) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = value, group = name,
                colour = name, fill = name
                )
            ) +
        ggplot2::geom_density(
            color = "white",
            alpha = 0.6,
            adjust = 5,
            show.legend = FALSE
            ) +
        ggplot2::geom_label(
            data = . %>% dplyr::summarise(m = unique(imag_rt_median) ),
            ggplot2::aes(x = m, y = 0, label = round(m, 3) ),
            position = ggplot2::position_nudge(y = 0.01),
            size = 4,
            inherit.aes = FALSE
            ) +
        ggplot2::geom_label(
            data = . %>% dplyr::summarise(m = unique(imag_mt_median) ),
            ggplot2::aes(x = m, y = 0, label = round(m, 3) ),
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

    # combining and returning the plots
    if (method == "functions") {

        p1

    } else if (method == "distributions") {

        p2

    } else if (method == "both") {

        p1 + p2 +
            patchwork::plot_layout(guides = "collect") &
            ggplot2::theme(legend.position = "bottom")

    }

}
