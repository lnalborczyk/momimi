#' Fitting the model
#'
#' Fitting the two versions of the model.
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
#' @param model_version Version of the model ("TMM" or "PIM").
#' @param method Optimisation method.
#' @param maxit Maximum number of iterations.
#'
#' @return The optimised parameter values and further convergence informatio
#'
#' @importFrom magrittr %>%
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
        nstudies = 200, initial_pop_constraints = FALSE,
        error_function, model_version,
        method = c(
            "SANN", "GenSA", "pso", "hydroPSO", "DEoptim",
            "Nelder-Mead", "BFGS", "L-BFGS-B", "bobyqa", "nlminb",
            "all_methods", "optimParallel"
            ),
        maxit = 1e2
        ) {

    # some tests for variable types
    stopifnot("data must be a dataframe..." = is.data.frame(data) )
    stopifnot("nsims must be a numeric..." = is.numeric(nsims) )
    stopifnot("nstudies must be a numeric..." = is.numeric(nstudies) )

    # method should be one of above
    method <- match.arg(method)

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

            # function for generating plausible starting values
            # source (file = paste0("utils/hypercube_sampling_while_", model_version, ".R") )

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
            fn = loss,
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
                trace = TRUE,
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
                packages = c("DEoptim", "tidyverse", "lhs", "momimi")
                # parVar = c("model", "loss_function")
                )
            )

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
