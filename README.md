
# Models of Motor Inhibition during Motor Imagery

[![GitHub repo
size](https://img.shields.io/github/repo-size/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)
[![GitHub last
update](https://img.shields.io/github/last-commit/lnalborczyk/momimi?color=brightgreen&logo=github)](https://github.com/lnalborczyk/momimi)

The `momimi` package implements the “threshold modulation model”
(TMM3/TMM4) and the “parallel inhibition model” (PIM) of motor
inhibition during motor imagery and provides fitting and plotting
utilities.

## Installation

You can install the latest stable version of the package from GitHub
with:

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
    amplitude_activ = 0.8, peak_time_activ = log(0.5), curvature_activ = 0.5,
    model_version = "TMM3",
    bw_noise = 0.05,
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
parameter values and then fit the model (below, the 3-parameter version
of the TMM) to these data to try recovering the original parameter
values.

``` r
# plausible "true" parameter values in the TMM3
# parameters are the motor execution threshold, the peak time, and the curvature
true_pars <- c(1.1, 0.5, 0.4)

# simulating data using these parameter values
simulated_data <- simulating(
    nsims = 200,
    nsamples = 3000,
    true_pars = true_pars,
    action_mode = "imagined",
    model_version = "TMM3"
    )

# displaying the first ten rows of these data
head(x = simulated_data, n = 10)
#>    reaction_time movement_time action_mode
#> 1      0.3467196     0.3468904    imagined
#> 2      0.3225490     0.4754676    imagined
#> 3      0.3149780     0.6991953    imagined
#> 4      0.2853013     0.7058372    imagined
#> 5      0.2950119     0.3818491    imagined
#> 6      0.3248526     0.4045011    imagined
#> 7      0.4005530     0.5672460    imagined
#> 8      0.3181484     0.4462748    imagined
#> 9      0.2860371     0.3530964    imagined
#> 10     0.3672168     0.4728014    imagined
```

We fit the model and use realistic constraints (e.g., the RT/MT should
be no less than 0.1s and no more than 2 seconds) on the initial
parameter values (by setting `initial_pop_constraints = TRUE`) to
facilitate convergence (fitting can take a while).

``` r
# fitting the model
results <- fitting(
    data = simulated_data,
    nsims = 200,
    error_function = "g2",
    method = "DEoptim",
    model_version = "TMM3",
    lower_bounds = c(1, 0.25, 0.1),
    upper_bounds = c(2, 1.25, 0.6),
    initial_pop_constraints = TRUE,
    maxit = 100
    )
```

``` r
# fitting summary
summary(results)
#> 
#> ***** summary of DEoptim object ***** 
#> best member   :  1.01658 0.50823 0.38563 
#> best value    :  0.01212 
#> after         :  100 generations 
#> fn evaluated  :  23129 times 
#> *************************************
```

We can then plot the underlying (latent) function(s). Note that this
returns a `ggplot2` object that can be subsequently modified.

``` r
plot(x = results, method = "latent", model_version = "TMM3") +
    labs(title = "Example of latent activation function")
```

<img src="man/figures/README-latent-1.png" width="75%" />

We can also do “predictive checks” by comparing the data used to fit the
model to data simulated from the model using the estimated parameter
values.

``` r
plot(
    x = results, original_data = simulated_data,
    method = "ppc", model_version = "TMM3", action_mode = "imagined"
    )
```

<img src="man/figures/README-ppc-1.png" width="75%" />

We can also do another form of “predictive check” by comparing the
quantile values of the original and simulated data.

``` r
plot(
    x = results, original_data = simulated_data,
    method = "quantiles", model_version = "TMM3", action_mode = "imagined"
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

Nalborczyk, L., Longcamp, M., Gajdos, T., Servant, M., & Alario, F.-X.
(under review). Towards formal models of inhibitory mechanisms involved
in motor imagery: A commentary on Bach, Frank, & Kunde (2022). Preprint
available at <https://psyarxiv.com/tz6x2/>.

## Getting help

If you encounter a bug or have a question please file an issue with a
minimal reproducible example on
[GitHub](https://github.com/lnalborczyk/momimi/issues).
