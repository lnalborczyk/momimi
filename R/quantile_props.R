#' Find quantile proportions
#'
#' Computes the proportion of observations within quantiles
#'
#' @param x Numeric, vector of observations (RTs or MTs)
#' @param quants Numeric, vector of quantiles bounds
#'
#' @return Vector of proportions per quantile.
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

quantiles_props <- function (x, quants) {

    quants2 <- as.numeric(c(0, quants, Inf) )
    quants_props <- as.numeric(table(cut(x, quants2) ) ) / length(x)

    return (quants_props)

}
