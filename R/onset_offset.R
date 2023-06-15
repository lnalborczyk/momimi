#' Finding onset and offset
#'
#' Finding onset and offset of the activation or balance functions (according to the model and some threshold).
#'
#' @param amplitude_activ Numeric, amplitude of the activation function.
#' @param amplitude_inhib Numeric, amplitude of the inhibition function.
#' @param peak_time_activ Numeric, peak time of the activation function.
#' @param peak_time_inhib Numeric, peak time of the inhibition function.
#' @param curvature_activ Numeric, curvature of the activation function.
#' @param curvature_inhib Numeric, curvature of the inhibition function.
#' @param thresh Numeric, a threshold defining onset and offset.
#' @param model_version Character, threshold modulation model ("TMM") or parallel inhibition model ("PIM").
#'
#' @return A vector containing the onset (i.e., the lowest value above the threshold) and the offset (i.e., the highest value above the threshold).
#'
#' @examples
#' \dontrun{
#' onset_offset(
#'     amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
#'     thresh = 1, model_version = "TMM"
#'     )
#' }
#'
#' @author Ladislas Nalborczyk \email{ladislas.nalborczyk@@gmail.com}.
#'
#' @export

onset_offset <- function (
        amplitude_activ = 1.5, peak_time_activ = 0, curvature_activ = 0.4,
        amplitude_inhib = 1.5, peak_time_inhib = 0, curvature_inhib = 0.6,
        thresh, model_version = c("TMM", "PIM")
        ) {

    # model_version should be one of above
    model_version <- match.arg(model_version)

    if (model_version == "TMM") {

        onset <- exp(peak_time_activ - sqrt(-2 * curvature_activ^2 * log(thresh / amplitude_activ) ) )
        offset <- exp(peak_time_activ + sqrt(-2 * curvature_activ^2 * log(thresh / amplitude_activ) ) )

        return (c(onset, offset) )

        } else if (model_version == "PIM") {

            a <- curvature_inhib^2 - curvature_activ^2
            b <- 2 * (curvature_activ^2 * peak_time_inhib - curvature_inhib^2 * peak_time_activ)
            c <- curvature_activ^2 * peak_time_inhib^2 - curvature_inhib^2 * peak_time_activ^2 - 2 * curvature_activ^2 * curvature_inhib^2 * (log(amplitude_activ / amplitude_inhib) - log(thresh) )
            onset <- exp((-b - sqrt(b^2 - 4 * a * c) ) / (2 * a) )
            offset <- exp((-b + sqrt(b^2 - 4 * a * c) ) / (2 * a) )

            return (c(onset, offset) )

        }

}
