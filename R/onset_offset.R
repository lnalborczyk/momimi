#' Finding onset and offset
#'
#' Finding onset and offset of the activation function (according to the model and some threshold).
#'
#' @param peak_time Numeric, peak time of the activation function.
#' @param curvature Numeric, curvature of the activation function.
#' @param thresh Numeric, a threshold defining onset and offset.
#'
#' @return A vector containing the onset (i.e., the lowest value above the threshold) and the offset (i.e., the highest value above the threshold).
#'
#' @examples
#' \dontrun{
#' onset_offset(peak_time = 0.5, curvature = 0.4, thresh = 0.5)
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

onset_offset <- function (peak_time = 0.5, curvature = 0.4, thresh = 0.5) {

    onset <- exp(peak_time - sqrt(-2 * curvature^2 * log(thresh) ) )
    offset <- exp(peak_time + sqrt(-2 * curvature^2 * log(thresh) ) )

    # direct formula for mt
    # mt <- 2 * exp(peak_time_activ) * sinh(curvature_activ * sqrt(2 * log(1 / thres) ) )

    return (c(onset, offset) )

}
