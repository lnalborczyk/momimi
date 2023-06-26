#' Fitting the model
#'
#' Fitting the "threshold modulation model" (TMM) and the "parallel inhibition model" (PIM) of motor inhibition during motor imagery.
#'
#' @param par Numeric, vector of parameter values.
#' @param data Dataframe, data to be used for fitting the model.
#' @param nsims Numeric, number of studies to be simulated.
#' @param par_names Character, vector of parameter names.
#' @param lower_bounds Numeric, vector of lower bounds for parameters.
#' @param upper_bounds Numeric, vector of upper bounds for parameters.
#' @param nstudies Numeric, number of starting values in the LHS.
#' @param initial_pop_constraints Boolean, whether to use additional constraints when sampling initial parameter values.
#' @param error_function Character, error function to be used when fitting the model.
#' @param model_version Character, version of the model ("TMM" or "PIM").
#' @param method Character, optimisation method (DEoptim seems to work best).
#' @param maxit Numeric, maximum number of iterations.
#' @param verbose Boolean, whether to print progress during fitting.
#'
#' @return The optimised parameter values and further convergence information.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' plausible "true" parameter values in the TMM
#' true_pars <- c(1.1, 0.5, 0.3, 1.25)
#'
#' # simulating data using these parameter values
#' simulated_data <- simulating(
#'     nsims = 200,
#'     nsamples = 2000,
#'     true_pars = true_pars,
#'     action_mode = "imagined",
#'     model_version = "TMM"
#'     )
#'
#' # fitting the model
#' results <- fitting(
#'     data = simulated_data,
#'     nsims = 200,
#'     error_function = "g2",
#'     method = "DEoptim",
#'     model_version = "TMM",
#'     lower_bounds = c(0, 0.25, 0.1, 1),
#'     upper_bounds = c(2, 1.25, 0.6, 2),
#'     initial_pop_constraints = TRUE,
#'     maxit = 100
#'     )
#'
#' # fitting summary
#' summary(results)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

# fitting the model
fitting <- function (
        par, data,
        nsims = NULL,
        par_names = NULL,
        lower_bounds, upper_bounds,
        nstudies = 200,
        initial_pop_constraints = FALSE,
        error_function = c("g2", "rmse", "sse", "wsse", "ks"),
        model_version = c("TMM", "PIM"),
        method = c(
            "SANN", "GenSA", "pso", "hydroPSO", "DEoptim",
            "Nelder-Mead", "BFGS", "L-BFGS-B", "bobyqa", "nlminb",
            "all_methods", "optimParallel"
            ),
        maxit = 100, verbose = TRUE
        ) {

    # defining parameter names according to the chosen model (if null)
    if (is.null(par_names) ) {

        if (model_version == "TMM") {

            par_names <- c("amplitude_activ", "peak_time_activ", "curvature_activ", "exec_threshold")

        } else if (model_version == "PIM") {

            par_names <- c("amplitude_ratio", "peak_time", "curvature_activ", "curvature_inhib")

        }

    }

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nstudies must be a numeric..." = is.numeric(nstudies) )

    # testing whether only 4 pars have been specified
    stopifnot("par_names must be a numeric of length 4..." = length(par_names) == 4)
    stopifnot("lower_bounds must be a numeric of length 4..." = length(lower_bounds) == 4)
    stopifnot("upper_bounds must be a numeric of length 4..." = length(upper_bounds) == 4)

    # method should be one of above
    method <- match.arg(method)

    # model_version should be one of above
    model_version <- match.arg(model_version)

    if (method == "SANN") {

        fit <- stats::optim(
            par = par,
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            method = method,
            control = list(maxit = maxit, trace = 2)
            )

    } else if (method == "GenSA") {

        fit <- GenSA::GenSA(
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            par = par,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = maxit, verbose = TRUE)
            )

    } else if (method == "pso") {

        fit <- pso::psoptim(
            fn = loss,
            data = data,
            par = par,
            nsims = nsims,
            error_function = error_function,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = maxit, trace = 2, trace.stats = TRUE)
            )

    } else if (method == "hydroPSO") {

        fit <- hydroPSO::hydroPSO(
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            par = par,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(
                maxit = maxit,
                verbose = TRUE,
                # using all available cores
                parallel = "parallel"
                )
            )

    } else if (method == "DEoptim") {

        if (initial_pop_constraints == TRUE) {

            # generating plausible starting values
            lhs_initial_pop <- generating_initialpop(
                nstudies = nstudies,
                action_mode = unique(data$action_mode),
                par_names = par_names,
                lower_bounds = lower_bounds, upper_bounds = upper_bounds,
                model_version = model_version
                )

        } else {

            # initialising an empty dataframe
            lhs_initial_pop_df <- data.frame(matrix(data = NA, nrow = nstudies, ncol = length(lower_bounds) ) )

            # populating it with hypercube samples
            for (i in 1:length(lower_bounds) ){

                lhs_initial_pop_df[, i] <- tgp::lhs(n = nstudies, rect = c(lower_bounds[i], upper_bounds[i]) )[, 1]

            }

            # converting to a matrix (as requested by deoptim)
            lhs_initial_pop <- as.matrix(lhs_initial_pop_df)

        }

        # starting the optimisation
        fit <- DEoptim::DEoptim(
            fn = momimi::loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            model_version = model_version,
            lower = lower_bounds,
            upper = upper_bounds,
            control = DEoptim::DEoptim.control(
                # maximum number of iterations
                itermax = maxit,
                # printing progress iteration
                trace = verbose,
                # printing progress every 10 iterations
                # trace = 10,
                # defines the differential evolution strategy (defaults to 2)
                # 1: DE / rand / 1 / bin (classical strategy)
                # 2: DE / local-to-best / 1 / bin (default)
                # 3: DE / best / 1 / bin with jitter
                # 4: DE / rand / 1 / bin with per-vector-dither
                # 5: DE / rand / 1 / bin with per-generation-dither
                # 6: DE / current-to-p-best / 1
                strategy = 3,
                # value to reach (defaults to -Inf)
                VTR = 0,
                # number of population members (by default 10*length(lower) )
                # NP = 200,
                NP = nrow(lhs_initial_pop),
                # F is the mutation constant (defaults to 0.8)
                F = 0.9,
                # crossover probability (recombination) (defaults to 0.5)
                CR = 0.95,
                # c controls the speed of the crossover adaptation
                # when strategy = 6 (defaults to 0)
                # c = 0.1,
                # proportion of best solutions to use in the mutation
                # when strategy = 6 (defaults to 0.2)
                # p = 0.1,
                # defining the initial population using lhs
                initialpop = lhs_initial_pop,
                # when to stop optimisation
                reltol = 1e-6,
                # number of iteration after which to stop the optimisation
                # if there is no improvement
                steptol = 1000,
                # using all available cores
                parallelType = "parallel",
                # packages = c("DEoptim", "tidyverse", "lhs", "momimi")
                packages = c("DEoptim", "tidyverse", "tgp", "momimi")
                )
            )

        # setting the class of the resulting object
        class(fit) <- c("DEoptim_momimi", "DEoptim", "data.frame")

    } else if (method %in% c("Nelder-Mead", "BFGS", "L-BFGS-B", "bobyqa", "nlminb") ) {

        fit <- optimx::optimx(
            par = par,
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            method = method,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = maxit, trace = 6)
            )

    } else if (method == "all_methods") {

        fit <- optimx::optimx(
            par = par,
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = maxit, trace = 2, all.methods = TRUE)
            )

    } else if (method == "optimParallel") {

        # using half of all available cores by default
        cl <- parallel::makeCluster(parallel::detectCores() / 2)

        # defining this as the default cluster
        parallel::setDefaultCluster(cl = cl)

        # loading the tidyverse package on each cluster
        parallel::clusterEvalQ(cl, library(tidyverse) )

        # loading model and loss function
        # clusterExport(cl, c("model", "loss_function") )

        # parallel optimisation
        fit <- optimParallel::optimParallel(
            par = par,
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            lower = lower_bounds,
            upper = upper_bounds,
            verbose = TRUE
            )

        # stopping the cluster
        parallel::stopCluster(cl)

    }

    return (fit)

}

#' @export

plot.DEoptim_momimi <- function (
        x, original_data,
        method = c("ppc", "latent"),
        action_mode = c("executed", "imagined"),
        model_version = c("TMM", "PIM"),
        ...
        ) {

    # some tests
    method <- match.arg(method)
    action_mode <- match.arg(action_mode)
    model_version <- match.arg(model_version)

    # retrieving estimated pars
    estimated_pars <- as.numeric(x$optim$bestmem)

    if (method == "ppc") {

        # simulating data using these parameter values
        simulating(
            nsims = 200,
            nsamples = 2000,
            true_pars = estimated_pars,
            action_mode = action_mode,
            model_version = model_version
            ) %>%
            # removing NAs or aberrant simulated data
            stats::na.omit() %>%
            dplyr::filter(.data$reaction_time < 3 & .data$movement_time < 3) %>%
            tidyr::pivot_longer(cols = .data$reaction_time:.data$movement_time) %>%
            ggplot2::ggplot(ggplot2::aes(x = .data$value, colour = .data$name, fill = .data$name) ) +
            ggplot2::geom_density(
                data = original_data %>% tidyr::pivot_longer(cols = .data$reaction_time:.data$movement_time),
                color = "white",
                position = "identity",
                alpha = 0.5, show.legend = FALSE
                ) +
            ggplot2::geom_density(size = 1, fill = NA, show.legend = FALSE) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                title = "Observed and simulated distributions of RTs/MTs",
                x = "Reaction/Movement time (in seconds)", y = "Probability density"
                )

    } else if (method == "latent") {

        if (model_version == "TMM") {

            par_names <- c("amplitude_activ", "peak_time_activ", "curvature_activ", "exec_threshold")

            parameters_estimates_summary <- paste(as.vector(rbind(
                paste0(par_names, ": "),
                paste0(as.character(round(estimated_pars, 3) ), "\n")
                ) ), collapse = "") %>% stringr::str_sub(end = -2)

            model(
                nsims = 500, nsamples = 2000,
                exec_threshold = estimated_pars[4] * estimated_pars[1],
                imag_threshold = 0.5 * estimated_pars[4] * estimated_pars[1],
                amplitude_activ = estimated_pars[1],
                peak_time_activ = log(estimated_pars[2]),
                curvature_activ = estimated_pars[3],
                model_version = "TMM",
                full_output = TRUE
                ) %>%
                tidyr::pivot_longer(cols = .data$activation) %>%
                ggplot2::ggplot(
                    ggplot2::aes(
                        x = .data$time, y = .data$value,
                        group = interaction(.data$sim, .data$name)
                        )
                    ) +
                ggplot2::geom_hline(
                    yintercept = estimated_pars[4] * estimated_pars[1],
                    linetype = 2
                    ) +
                ggplot2::geom_hline(
                    yintercept = 0.5 * estimated_pars[4] * estimated_pars[1],
                    linetype = 2
                    ) +
                # plotting average
                ggplot2::stat_summary(
                    # ggplot2::aes(group = .data$name, colour = .data$name),
                    ggplot2::aes(group = .data$name),
                    fun = "median", geom = "line",
                    linewidth = 1, alpha = 1,
                    show.legend = FALSE
                    ) +
                # displaying estimated parameter values
                ggplot2::annotate(
                    geom = "label",
                    x = Inf, y = Inf,
                    hjust = 1, vjust = 1,
                    label = parameters_estimates_summary,
                    family = "Courier"
                    ) +
                ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
                ggplot2::labs(
                    title = "Latent functions",
                    x = "Time within a trial (in seconds)",
                    y = "Activation/inhibition (a.u.)",
                    colour = "",
                    fill = ""
                    )

        } else if (model_version == "PIM") {

            par_names <- c("amplitude_ratio", "peak_time", "curvature_activ", "curvature_inhib")

            parameters_estimates_summary <- paste(as.vector(rbind(
                paste0(par_names, ": "),
                paste0(as.character(round(estimated_pars, 3) ), "\n")
                ) ), collapse = "") %>% stringr::str_sub(end = -2)

            model(
                nsims = 500, nsamples = 2000,
                exec_threshold = 1, imag_threshold = 0.5,
                amplitude_activ = 1.5,
                peak_time_activ = log(estimated_pars[2]),
                curvature_activ = estimated_pars[3],
                amplitude_inhib = 1.5 / estimated_pars[1],
                peak_time_inhib = log(estimated_pars[2]),
                curvature_inhib = estimated_pars[4] * estimated_pars[3],
                model_version = "PIM",
                full_output = TRUE
                ) %>%
                tidyr::pivot_longer(cols = .data$activation:.data$balance) %>%
                ggplot2::ggplot(
                    ggplot2::aes(
                        x = .data$time, y = .data$value,
                        group = interaction(.data$sim, .data$name),
                        colour = .data$name
                        )
                    ) +
                ggplot2::geom_hline(yintercept = 1, linetype = 2) +
                ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
                # plotting average
                ggplot2::stat_summary(
                    ggplot2::aes(group = .data$name, colour = .data$name),
                    fun = "median", geom = "line",
                    linewidth = 1, alpha = 1,
                    show.legend = TRUE
                    ) +
                # displaying estimated parameter values
                ggplot2::annotate(
                    geom = "label",
                    x = Inf, y = Inf,
                    hjust = 1, vjust = 1,
                    label = parameters_estimates_summary,
                    family = "Courier"
                    ) +
                ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
                ggplot2::labs(
                    title = "Latent activation, inhibition, and balance functions",
                    x = "Time within a trial (in seconds)",
                    y = "Activation/inhibition (a.u.)",
                    colour = "",
                    fill = ""
                    )

        }

    }

}
