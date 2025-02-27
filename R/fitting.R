#' Fitting the model
#'
#' Fitting the "threshold modulation model" (TMM) of motor inhibition during motor imagery.
#'
#' @param par Numeric, vector of parameter values.
#' @param data Dataframe, data to be used for fitting the model.
#' @param nsims Numeric, number of studies to be simulated.
#' @param par_names Character, vector of parameter names.
#' @param lower_bounds Numeric, vector of lower bounds for parameters.
#' @param upper_bounds Numeric, vector of upper bounds for parameters.
#' @param nstudies Numeric, number of starting values in the LHS.
#' @param initial_pop_constraints Boolean, whether to use additional constraints when sampling initial parameter values.
#' @param rt_contraints Numeric, vector of length 2 specifying the min and max RT (in seconds).
#' @param mt_contraints Numeric, vector of length 2 specifying the min and max MT (in seconds).
#' @param error_function Character, error function to be used when fitting the model.
#' @param uncertainty Numeric, indicates how noise is introduced in the system.
#' @param method Character, optimisation method (DEoptim seems to work best).
#' @param cluster Character, existing parallel cluster object. If provided, overrides + specified parallelType.
#' @param grid_resolution Numeric, resolution of the grid when method = "grid_search".
#' @param maxit Numeric, maximum number of iterations.
#' @param verbose Boolean, whether to print progress during fitting.
#'
#' @return The optimised parameter values and further convergence information.
#'
#' @importFrom magrittr %>%
#' @importFrom utils head
#' @importFrom stats density
#'
#' @examples
#' \dontrun{
#' plausible "true" parameter values in the TMM (with 4 free parameters)
#' true_pars <- c(1.1, 0.5, 0.3, 0.08)
#'
#' # simulating data using these parameter values
#' simulated_data <- simulating(
#'     nsims = 200,
#'     true_pars = true_pars,
#'     action_mode = "imagined"
#'     )
#'
#' # fitting the model (to keep a parameter constant, assign the same value
#' # to both upper and lower...)
#' results <- fitting(
#'     data = simulated_data,
#'     nsims = 200,
#'     error_function = "g2",
#'     method = "DEoptim",
#'     lower_bounds = c(1, 0.25, 0.3, 0.01),
#'     upper_bounds = c(2, 1.25, 0.3, 0.20),
#'     initial_pop_constraints = TRUE,
#'     maxit = 20
#'     )
#'
#' # fitting summary
#' summary(results)
#'
#' # plotting the results
#' plot(results, method = "latent")
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
        rt_contraints = c(0.1, 2), mt_contraints = c(0.1, 2),
        error_function = c("g2", "rmse", "sse", "wsse", "ks"),
        uncertainty = c("par_level", "func_level", "diffusion"),
        method = c(
            "SANN", "GenSA", "pso", "DEoptim",
            "Nelder-Mead", "BFGS", "L-BFGS-B", "bobyqa", "nlminb",
            "all_methods", "optimParallel", "grid_search"
            ),
        cluster = NULL,
        grid_resolution = 0.01,
        maxit = 100, verbose = TRUE
        ) {

    # defining parameter names according to the chosen model (if null)
    if (is.null(par_names) ) {

        par_names <- c("exec_threshold", "peak_time", "curvature", "bw_noise")

    }

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nstudies must be a numeric..." = is.numeric(nstudies) )

    # testing whether 4 pars have been specified
    stopifnot("par_names must be a numeric of length 4..." = length(par_names) == 4)
    stopifnot("lower_bounds must be a numeric of length 4..." = length(lower_bounds) == 4)
    stopifnot("upper_bounds must be a numeric of length 4..." = length(upper_bounds) == 4)

    # method should be one of above
    method <- match.arg(method)

    # uncertainty should be one of above
    uncertainty <- match.arg(uncertainty)

    if (method == "SANN") {

        fit <- stats::optim(
            par = par,
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            uncertainty = uncertainty,
            method = method,
            control = list(maxit = maxit, trace = 2)
            )

    } else if (method == "GenSA") {

        fit <- GenSA::GenSA(
            fn = loss,
            data = data,
            nsims = nsims,
            error_function = error_function,
            uncertainty = uncertainty,
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
            uncertainty = uncertainty,
            lower = lower_bounds,
            upper = upper_bounds,
            control = list(maxit = maxit, trace = 2, trace.stats = TRUE)
            )

    } else if (method == "DEoptim") {

        if (initial_pop_constraints == TRUE) {

            # generating plausible starting values
            lhs_initial_pop <- generating_initialpop(
                nstudies = nstudies,
                action_mode = unique(data$action_mode),
                par_names = par_names,
                lower_bounds = lower_bounds,
                upper_bounds = upper_bounds,
                uncertainty = uncertainty,
                rt_contraints = rt_contraints, mt_contraints = mt_contraints
                )

            # printing the first few starting values (sanity check)
            head(x = lhs_initial_pop)

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
            uncertainty = uncertainty,
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
                CR = 0.9,
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
                # defining the package to be imported on each parallel core
                packages = c("DEoptim", "dplyr", "tidyr", "tgp", "momimi"),
                # defining the cluster
                cluster = cluster
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
            uncertainty = uncertainty,
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
            uncertainty = uncertainty,
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
            uncertainty = uncertainty,
            lower = lower_bounds,
            upper = upper_bounds,
            verbose = TRUE
            )

        # stopping the cluster
        parallel::stopCluster(cl)

    } else if (method == "grid_search") {

        # using all available cores
        future::plan(future::multisession(workers = parallel::detectCores() ) )

        # defining the grid of parameter values
        parameters_grid <- tidyr::crossing(
            a = seq(from = lower_bounds[1], to = upper_bounds[1], by = grid_resolution),
            b = seq(from = lower_bounds[2], to = upper_bounds[2], by = grid_resolution),
            c = seq(from = lower_bounds[3], to = upper_bounds[3], by = grid_resolution),
            d = seq(from = lower_bounds[4], to = upper_bounds[4], by = grid_resolution)
            )

        # warning the user about the number of simulation to evaluate...
        message(
            paste(
                "momimi will now explore",
                nrow(parameters_grid),
                "combinations of parameters values, so please adjust your expectations accordingly..."
                )
            )

        # setting up the progress bar (cf. https://progressr.futureverse.org)
        progressr::handlers(global = TRUE)

        # initialising the progress bar
        p <- progressr::progressor(steps = nrow(parameters_grid) )

        # computing the error for many possible parameters values
        error_surface <- parameters_grid %>%
            dplyr::mutate(
                error = future.apply::future_apply(
                    X = ., MARGIN = 1,
                    FUN = function (x, ...) {
                        p(sprintf("x=%g", x) )
                        momimi::loss(x, data = data, nsims = nsims)
                        },
                    # FUN = momimi::loss,
                    # data = data,
                    # nsims = nrow(data),
                    future.seed = NULL
                    )
                )

        # explicitly closing multisession workers by switching plan
        future::plan(future::sequential)

        # returning the error surface
        fit <- error_surface

    }

    return (fit)

}

#' @export

plot.DEoptim_momimi <- function (
        x, original_data,
        method = c("ppc", "quantiles", "latent", "optimisation"),
        action_mode = c("executed", "imagined"),
        uncertainty = c("par_level", "func_level", "diffusion"),
        nsims = 500, nsamples = 3000,
        max_rt = 5, max_mt = 5,
        ...
        ) {

    # matching arguments
    method <- match.arg(method)
    action_mode <- match.arg(action_mode)
    uncertainty <- match.arg(uncertainty)

    # retrieving estimated pars
    estimated_pars <- as.numeric(x$optim$bestmem)

    if (method == "ppc") {

        # simulating data using these parameter values
        simulated_data <- simulating(
            nsims = nsims,
            nsamples = nsamples,
            true_pars = estimated_pars,
            action_mode = action_mode,
            uncertainty = uncertainty
            ) %>%
            # removing NAs or aberrant simulated data
            stats::na.omit() %>%
            dplyr::filter(.data$reaction_time < max_rt & .data$movement_time < max_mt) %>%
            tidyr::pivot_longer(cols = .data$reaction_time:.data$movement_time)

        # computing the "optimal" bandwidth
        optimal_bw <- simulated_data %>%
            dplyr::pull(.data$value) %>%
            stats::bw.nrd()


        # sanity check
        # cat(optimal_bw)

        # plotting
        simulated_data %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = .data$value,
                    colour = .data$name,
                    fill = .data$name
                    )
                ) +
            ggplot2::geom_histogram(
                data = original_data %>%
                    tidyr::pivot_longer(
                        cols = .data$reaction_time:.data$movement_time
                        ),
                ggplot2::aes(y = ggplot2::after_stat(density), fill = .data$name),
                position = "identity",
                alpha = 0.6,
                colour = "white",
                binwidth = optimal_bw,
                show.legend = FALSE
                ) +
            ggplot2::geom_density(size = 1, fill = NA, show.legend = FALSE) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                title = "Observed and simulated distributions of RTs/MTs",
                x = "Reaction/Movement time (s)", y = "Probability density"
                )

    } else if (method == "quantiles") {

        # simulating data using these parameter values
        simulated_data <- simulating(
            nsims = nsims,
            nsamples = nsamples,
            true_pars = estimated_pars,
            action_mode = action_mode,
            uncertainty = uncertainty
            ) %>%
            # removing NAs or aberrant simulated data
            stats::na.omit() %>%
            dplyr::filter(.data$reaction_time < max_rt & .data$movement_time < max_mt)

        # what quantiles should we look at?
        quantile_probs <- seq(0.1, 0.9, 0.1)

        # binding original and simulated data and plotting it
        dplyr::bind_rows(
            original_data %>% dplyr::mutate(type = "observed"),
            simulated_data %>% dplyr::mutate(type = "simulated")
            ) %>%
            tidyr::pivot_longer(names_to = "measure", cols = .data$reaction_time:.data$movement_time) %>%
            dplyr::group_by(.data$type, .data$measure) %>%
            dplyr::reframe(tibble::enframe(stats::quantile(x = .data$value, probs = quantile_probs, na.rm = TRUE, names = TRUE) ) ) %>%
            dplyr::ungroup() %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = .data$name, y = .data$value,
                    group = interaction(.data$type, .data$measure),
                    colour = .data$measure, fill = .data$measure,
                    shape = .data$type
                    )
                ) +
            ggplot2::geom_line(alpha = 0.5, linetype = 3) +
            ggplot2::geom_point(
                size = 4,
                alpha = 0.9,
                show.legend = TRUE
                ) +
            ggplot2::theme_bw(base_size = 12, base_family = "Open Sans") +
            ggplot2::labs(
                x = "Percentile", y = "Time (s)",
                colour = "", fill = "", shape = ""
                )

    } else if (method == "latent") {

        par_names <- c("exec_threshold", "peak_time", "curvature", "bw_noise")

        parameters_estimates_summary <- paste(as.vector(rbind(
            paste0(par_names, ": "),
            paste0(as.character(round(estimated_pars, 3) ), "\n")
            ) ), collapse = "") %>% stringr::str_sub(end = -2)

        model(
            nsims = nsims,
            nsamples = nsamples,
            exec_threshold = estimated_pars[1],
            imag_threshold = 0.5 * estimated_pars[1],
            peak_time = log(estimated_pars[2]),
            curvature = estimated_pars[3],
            bw_noise = estimated_pars[4],
            uncertainty = uncertainty,
            full_output = TRUE
            ) %>%
            tidyr::pivot_longer(cols = .data$activation) %>%
            ggplot2::ggplot(
                ggplot2::aes(
                    x = .data$time, y = .data$value,
                    group = interaction(.data$sim, .data$name)
                    )
                ) +
            # plotting the motor execution and motor imagery thresholds
            ggplot2::geom_hline(
                yintercept = estimated_pars[1],
                linetype = 2
                ) +
            ggplot2::geom_hline(
                yintercept = 0.5 * estimated_pars[1],
                linetype = 2
                ) +
            # plotting average
            ggplot2::stat_summary(
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
                title = "Latent activation function",
                x = "Time within a trial (s)",
                y = "Activation (a.u.)",
                colour = "",
                fill = ""
                )

    } else if (method == "optimisation") {

        # plotting optimisation paths in parameter space
        optimisation_results <- data.frame(x$member$bestmemit) %>%
            dplyr::mutate(iteration = 1:nrow(.) )

        # dynamic 3D plot
        plotly::plot_ly(
            data = dplyr::distinct(optimisation_results),
            x = ~.data$par1, y = ~.data$par2, z = ~.data$par3
            ) %>%
            plotly::layout(
                scene = list(
                    xaxis = list(title = "Parameter 1"),
                    yaxis = list(title = "Parameter 2"),
                    zaxis = list(title = "Parameter 3")
                    )
                ) %>%
            plotly::add_trace(type = "scatter3d", mode = "markers+lines", color = ~iteration)

    }

}
