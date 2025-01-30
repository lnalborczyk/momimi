# Models of Motor Inhibition during Motor Imagery

[![GitHub repo
size](https://img.shields.io/github/repo-size/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub last
update](https://img.shields.io/github/last-commit/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)

The `momimi` package implements the “threshold modulation model” (TMM)
of motor imagery and provides fitting and plotting utilities in `R`.

## Installation

You can install the latest stable version of the package from GitHub
with:

``` r
# if needed, install.packages("devtools")
devtools::install_github(repo = "lnalborczyk/momimi", build_vignettes = TRUE)
```

## Usage

### Simulating and plotting data

We start by simulating some data (here, 100 observations or RTs and MTs)
given some values for the model’s parameters.

``` r
library(tidyverse)
library(momimi)

simulated_data <- model(
    nsims = 100, nsamples = 2000,
    exec_threshold = 1.1, imag_threshold = 0.55,
    peak_time = log(0.5), curvature = 0.5,
    full_output = TRUE
    )
```

We can plot the underlying activation function and the simulated
distributions of RTs and MTs.

``` r
# plotting the latent activation function
plot(x = simulated_data, method = "functions")
```

<img src="man/figures/README-plotting-1.png" width="75%" />

``` r

# plotting the distributions of RTs/MTs
plot(x = simulated_data, method = "distributions")
```

<img src="man/figures/README-plotting-2.png" width="75%" />

### Fitting the models

We can also use the model to generate realistic data from known
parameter values and then fit the model to these data to try recovering
the original parameter values.

``` r
# defining plausible "true" parameter values:
# the motor execution threshold, peak time, curvature, and between-trial variability (SD)
true_pars <- c(1.1, 0.5, 0.4, 0.1)

# simulating data using these parameter values
simulated_data <- simulating(
    nsims = 200,
    nsamples = 3000,
    true_pars = true_pars,
    action_mode = "imagined"
    )

# displaying the first ten rows of these data
head(x = simulated_data, n = 10)
#>    reaction_time movement_time action_mode
#> 1      0.2750805     0.4806944    imagined
#> 2      0.2963450     0.4955550    imagined
#> 3      0.2690255     0.5035881    imagined
#> 4      0.3284638     0.6321461    imagined
#> 5      0.3523583     0.3937026    imagined
#> 6      0.2613692     0.3791857    imagined
#> 7      0.3277495     0.6520622    imagined
#> 8      0.3158585     0.4149851    imagined
#> 9      0.3383590     0.3696656    imagined
#> 10     0.3368878     0.2612512    imagined
```

We fit the model and use realistic constraints (e.g., the RT/MT should
be no less than 0.1s and no more than 2 seconds) on the initial
parameter values (by setting `initial_pop_constraints = TRUE`) to
facilitate convergence (fitting can take a while).

``` r
# fitting the model
results <- fitting(
    # the data used to fit the model
    data = simulated_data,
    # the number of simulated trials in the model
    nsims = 200,
    # the g2 cost function is also used with the DDM
    error_function = "g2",
    # DEoptim seems the most efficient approach
    method = "DEoptim",
    # lower and upper bounds for the parameters values
    # NB: one can fix (i.e., not estimate) a parameter 
    # by setting the lower and upper bounds to the same value
    lower_bounds = c(1.0, 0.3, 0.4, 0.05),
    upper_bounds = c(1.4, 0.7, 0.4, 0.20),
    # should we generate sensible starting values?
    initial_pop_constraints = TRUE,
    nstudies = 200,
    # which constraints?
    rt_contraints = range(simulated_data$reaction_time),
    mt_contraints = range(simulated_data$movement_time),
    # maximum number of iteration when fitting the model (to be increased)
    maxit = 10
    )
```

``` r
# fitting summary
summary(results)
#> 
#> ***** summary of DEoptim object ***** 
#> best member   :  1.15826 0.49017 0.4 0.10122 
#> best value    :  0.02885 
#> after         :  10 generations 
#> fn evaluated  :  3190 times 
#> *************************************
```

We can then plot the underlying (latent) activation function. Note that
this returns a `ggplot2` object that can be subsequently modified.

``` r
plot(x = results, method = "latent") +
    labs(title = "Estimated latent activation function")
```

<img src="man/figures/README-latent-1.png" width="75%" />

We can also do “predictive checks” by comparing the data used to fit the
model to data simulated from the model using the estimated parameter
values.

``` r
plot(
    x = results, original_data = simulated_data,
    method = "ppc", action_mode = "imagined"
    )
```

<img src="man/figures/README-ppc-1.png" width="75%" />

We can also do another form of “predictive check” by comparing the
quantile values of the original and simulated data.

``` r
plot(
    x = results, original_data = simulated_data,
    method = "quantiles", action_mode = "imagined"
    )
```

<img src="man/figures/README-quantiles-1.png" width="100%" />

We can also visualise the trajectory in parameter space during
optimisation interactively using `plotly` (not shown below, but you
should try this).

``` r
plot(x = results, method = "optimisation")
```

## References

Nalborczyk, L., Longcamp, M., Gajdos, T., Servant, M. & Alario, F.‐X.
(in preparation). Modelling the onset and duration of imagined actions:
Assessing a novel algorithmic model of motor imagery.

Nalborczyk, L., Longcamp, M., Gajdos, T., Servant, M. & Alario, F.‐X.
(2024). Towards formal models of inhibitory mechanisms involved in motor
imagery: A commentary on Bach et al. (2022). Psychological Research,
1‐4. <https://doi.org/10.1007/s00426-023-01915-8>. Preprint available at
<https://psyarxiv.com/tz6x2/>.

## Getting help

If you encounter a bug or have a question please file an issue with a
minimal reproducible example on
[GitHub](https://github.com/lnalborczyk/momimi/issues).
