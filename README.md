
# Models of Motor Inhibition during Motor Imagery

[![GitHub repo
size](https://img.shields.io/github/repo-size/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub last
update](https://img.shields.io/github/last-commit/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub
downloads](https://img.shields.io/github/downloads/lnalborczyk/momimi/total?logo=github)](https://github.com/lnalborczyk/momimi)

The `momimi` package implements the “threshold modulation model”
(TMM3/TMM4) and the “parallel inhibition model” (PIM) of motor
inhibition during motor imagery and provides fitting and plotting
utilities.

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
    model_version = "TMM3",
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
# plausible "true" parameter values in the TMM
true_pars <- c(1.1, 0.5, 0.4)

# simulating data using these parameter values
simulated_data <- simulating(
    nsims = 200,
    nsamples = 2000,
    true_pars = true_pars,
    action_mode = "imagined",
    model_version = "TMM3"
    )

# displaying the first ten rows of these data
head(x = simulated_data, n = 10)
#>    reaction_time movement_time action_mode
#> 1      0.3191180     0.4563459    imagined
#> 2      0.3281299     0.4424445    imagined
#> 3      0.3166157     0.4672474    imagined
#> 4      0.3260946     0.4475316    imagined
#> 5      0.3178937     0.4536711    imagined
#> 6      0.3223118     0.4705083    imagined
#> 7      0.3159033     0.4603504    imagined
#> 8      0.3248285     0.4388826    imagined
#> 9      0.3280234     0.4430937    imagined
#> 10     0.3226074     0.4283027    imagined
```

We fit the model and use extra constraints on the initial parameter
values to facilitate convergence (fitting can take a while).

``` r
# fitting the model
results <- fitting(
    data = simulated_data,
    nsims = 500,
    error_function = "g2",
    method = "DEoptim",
    model_version = "TMM3",
    lower_bounds = c(1, 0.25, 0.1),
    upper_bounds = c(2, 1.25, 0.6),
    initial_pop_constraints = FALSE,
    maxit = 100
    )
```

``` r
# fitting summary
summary(results)
#> 
#> ***** summary of DEoptim object ***** 
#> best member   :  1.07338 0.4994 0.39119 
#> best value    :  0.00408 
#> after         :  100 generations 
#> fn evaluated  :  20200 times 
#> *************************************
```

We can then plot the underlying (latent) function(s).

``` r
plot(x = results, method = "latent", model_version = "TMM3") +
    labs(title = "Example of latent activation function")
```

<img src="man/figures/README-latent-1.png" width="75%" />

We can also do predictive checks by comparing the data used to fit the
model to data simulated from the model using the estimated parameter
values.

``` r
plot(
    x = results, original_data = simulated_data,
    method = "ppc", model_version = "TMM3", action_mode = "imagined"
    )
```

<img src="man/figures/README-ppc-1.png" width="75%" />

## References

Nalborczyk, L., Longcamp, M., Gajdos, T., Servant, M., & Alario, F.-X.
(*submitted*). Towards formal models of inhibitory mechanisms involved
in motor imagery: A commentary on Bach, Frank, & Kunde (2022). Preprint
available at <https://psyarxiv.com/tz6x2/>.

## Getting help

If you encounter a bug or have a question please file an issue with a
minimal reproducible example on
[GitHub](https://github.com/lnalborczyk/momimi/issues).
