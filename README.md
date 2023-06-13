
# Models of Motor Inhibition during Motor Imagery

[![Build
Status](https://travis-ci.org/lnalborczyk/momimi.svg?branch=master)](https://travis-ci.org/lnalborczyk/momimi)
[![GitHub repo
size](https://img.shields.io/github/repo-size/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub last
update](https://img.shields.io/github/last-commit/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub
downloads](https://img.shields.io/github/downloads/lnalborczyk/momimi/total&logo=github)](https://github.com/lnalborczyk/momimi)

The `momimi` package implements the “threshold modulation model” (TMM)
and the “parallel inhibition model” (PIM) of motor inhibition during
motor imagery and provides several fitting and plotting utilities.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github(repo = "lnalborczyk/momimi", build_vignettes = TRUE)
```

## Usage

### Simulating and plotting data

We start by simulating some data (here, 100 observations or RTs and
MTs).

``` r
library(tidyverse)
library(momimi)

simulated_data <- model(
    nsims = 100, nsamples = 2000,
    exec_threshold = 1, imag_threshold = 0.5,
    amplitude_activ = 0.8, peak_time_activ = log(0.5), curvature_activ = 0.4,
    model_version = "TMM",
    full_output = TRUE
    )
```

We can plot the underlying activation function and the implied
distributions of RTs and MTs.

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

### Fitting the models

We can also use the model to generate realistic data from known
parameter values and then fit the model to these data to try recovering
the original parameter values.

``` r
# plausible "true" parameter values
true_pars <- c(1.1, 0.5, 0.3, 1.25)

# simulating data using these parameter values
simulated_data <- model(
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
head(x = simulated_data, n = 10)
#>    reaction_time movement_time action_mode
#> 1      0.3712421     0.2942656    imagined
#> 2      0.3806787     0.2684245    imagined
#> 3      0.3715586     0.3110980    imagined
#> 4      0.3671511     0.3004650    imagined
#> 5      0.3724952     0.2979726    imagined
#> 6      0.3663551     0.3092010    imagined
#> 7      0.3721625     0.3078148    imagined
#> 8      0.3818147     0.2841863    imagined
#> 9      0.3673456     0.2945165    imagined
#> 10     0.3736252     0.2936363    imagined
```

``` r
# fitting the model
results <- fitting(
    data = simulated_data,
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
#> Iteration: 1 bestvalit: 2.484282 bestmemit:    0.397222    0.508942    0.250493    0.893974
#> Iteration: 2 bestvalit: 1.932035 bestmemit:    1.657351    0.512966    0.156453    0.349366
#> Iteration: 3 bestvalit: 1.932035 bestmemit:    1.657351    0.512966    0.156453    0.349366
#> Iteration: 4 bestvalit: 1.932035 bestmemit:    1.657351    0.512966    0.156453    0.349366
#> Iteration: 5 bestvalit: 1.470281 bestmemit:    1.975907    0.503800    0.149009    0.340395
#> Iteration: 6 bestvalit: 1.470281 bestmemit:    1.975907    0.503800    0.149009    0.340395
#> Iteration: 7 bestvalit: 1.470281 bestmemit:    1.975907    0.503800    0.149009    0.340395
#> Iteration: 8 bestvalit: 0.900814 bestmemit:    0.386890    0.509362    0.162787    0.352842
#> Iteration: 9 bestvalit: 0.900814 bestmemit:    0.386890    0.509362    0.162787    0.352842
#> Iteration: 10 bestvalit: 0.900814 bestmemit:    0.386890    0.509362    0.162787    0.352842
#> Iteration: 11 bestvalit: 0.900814 bestmemit:    0.386890    0.509362    0.162787    0.352842
#> Iteration: 12 bestvalit: 0.900814 bestmemit:    0.386890    0.509362    0.162787    0.352842
#> Iteration: 13 bestvalit: 0.617094 bestmemit:    0.987943    0.503503    0.175923    0.510871
#> Iteration: 14 bestvalit: 0.328873 bestmemit:    0.653311    0.503199    0.244202    0.990749
#> Iteration: 15 bestvalit: 0.328873 bestmemit:    0.653311    0.503199    0.244202    0.990749
#> Iteration: 16 bestvalit: 0.328873 bestmemit:    0.653311    0.503199    0.244202    0.990749
#> Iteration: 17 bestvalit: 0.328873 bestmemit:    0.653311    0.503199    0.244202    0.990749
#> Iteration: 18 bestvalit: 0.328873 bestmemit:    0.653311    0.503199    0.244202    0.990749
#> Iteration: 19 bestvalit: 0.270176 bestmemit:    0.584146    0.503199    0.244202    0.990749
#> Iteration: 20 bestvalit: 0.270176 bestmemit:    0.584146    0.503199    0.244202    0.990749

# fitting summary
summary(results)
#> 
#> ***** summary of DEoptim object ***** 
#> best member   :  0.58415 0.5032 0.2442 0.99075 
#> best value    :  0.27018 
#> after         :  20 generations 
#> fn evaluated  :  4200 times 
#> *************************************
```

## References

Nalborczyk, L., Longcamp, M., Gajdos, T., Servant, M., & Alario, F.-X.
(*to be submitted*). Towards formal models of inhibitory mechanisms
involved in motor imagery: A commentary on Bach, Frank, & Kunde (2022).

## Getting help

If you encounter a bug or have a question please file an issue with a
minimal reproducible example on
[GitHub](https://github.com/lnalborczyk/momimi/issues).
