
# Models of Motor Inhibition during Motor Imagery (momimi)

<!-- badges: start -->
<!-- badges: end -->

The goal of `momimi` is to provide utilities for fitting two models of
motor inhibition during motor imagery.

## Installation

You can install the development version of `momimi` from GitHub with:

``` r
remotes::install_github(
    repo = "https://github.com/lnalborczyk/momimi",
    dependencies = TRUE
    )
```

## Usage

We start by simulating some data (here, 100 observations or RTs and
MTs).

``` r
library(momimi)

simulated_data <- model(
    nsims = 100, nsamples = 2000,
    exec_threshold = 1, imag_threshold = 0.5,
    amplitude_activ = 0.8, peak_time_activ = log(0.5), curvature_activ = 0.4,
    model_version = "TMM",
    full_output = TRUE
    )
```

We can plot the data…

``` r
# plotting only the latent function(s)
plot(x = simulated_data, method = "functions")
```

<img src="man/figures/README-plotting-1.png" width="75%" />

``` r

# plotting only the distributions of RTs/MTs distributions
plot(x = simulated_data, method = "distributions")
```

<img src="man/figures/README-plotting-2.png" width="75%" />

We can also use the model to generate realistic data from known
parameter values and then fit the model to these data to try recover the
original parameter values.

``` r
# plausible "true" parameter values
true_pars <- c(1.1, 0.5, 0.3, 1.25)

# simulating data using these parameter values
df <- model(
    nsims = 200, nsamples = 2000,
    exec_threshold = true_pars[4] * true_pars[1],
    imag_threshold = 0.5 * true_pars[4] * true_pars[1],
    amplitude_activ = true_pars[1],
    peak_time_activ = log(true_pars[2]),
    curvature_activ = true_pars[3],
    model_version = "TMM",
    full_output = FALSE
    ) %>%
    mutate(action_mode = "imagined") %>%
    # keeping only the relevant columns
    dplyr::select(
        sim,
        reaction_time = paste0("onset_", substr(unique(.$action_mode), 1, 4) ),
        movement_time = paste0("mt_", substr(unique(.$action_mode), 1, 4) ),
        action_mode
        ) %>%
    distinct() %>%
    dplyr::select(-sim)

# displaying the first ten rows of these data
head(x = df, n = 10)
```

``` r
# fitting the model
results <- fitting(
    data = df,
    nsims = 200,
    error_function = "g2",
    method = "DEoptim",
    model_version = "TMM",
    par_names = c("amplitude_activ", "peak_time_activ", "curvature_activ", "exec_threshold"),
    lower_bounds = c(0, 0.5, 0, 0),
    upper_bounds = c(2, 1.5, 1, 1),
    nstudies = 200,
    initial_pop_constraints = FALSE,
    maxit = 20
    )

# fitting summary
summary(results)
```
