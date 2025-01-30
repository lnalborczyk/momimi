# momimi 0.0.1

* Adding a `NEWS.md` file to track changes.

# momimi 0.0.2

* Adding the 3-parameter version of TMM (free parameters are the relative motor execution threshold, the peak time, and the curvature).

# momimi 0.0.3

* Adding the possibility to include within-trial (diffusive) noise.
* Fixing the 3-parameter version of the TMM (TMM3).

# momimi 0.0.4

* Fixing the plotting method of `momimi_sim`.
* It is now possible to estimate the amount of between-trial variability (noise) in the TMM3 (which thus now has 4 free parameters).

# momimi 0.0.5

* Fixing some erroneous labelling in the plotting utilities.
* Adding the quantiles option in plotting methods for fitted objects (qq-plot).
* Adding the possibility to fit models with function-level noise or diffusive noise (in addition to only parameter-level noise). These two options are slower because it requires numerically finding the RT and MT.
* The 3-parameter version of the TMM (TMM3) now really has only 3 parameters, whereas the 4-parameter version now has the motor execution threshold, the peak time, the curvature, and the amount of between-trial variability as free parameters.
* It is now possible to fit the models with "brute force" (i.e., defining a discrete grid and looking for the minimum).

# momimi 0.0.6

* Fixing some erroneous labelling in the plotting utilities.
* Now returning the full error surface in `grid_search` fitting method.
* Removing constraints during fitting (i.e., exec_threshold was previously expected between 0.25 and 4).

# momimi 0.0.7

* Removing the parallel inhibition model (non-identifiable).
* Simplifying the threshold modulation model by dropping the amplitude parameter (fixed to 1).

# momimi 0.0.8

* Adding the possibility to specify the number of available cores in `fitting()`.

# momimi 0.0.9

* Fixing some typos.
* Fixing an error related to the presence of a `balance` variable (defined in a previous version of the model).
* Removed the 3-parameter version of the TMM: now the number of free parameters can simply be adjusted by modifying the `lower_bounds` and `upper_bounds` arguments when using the `fitting()` function.
