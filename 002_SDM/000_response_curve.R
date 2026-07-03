#' Generate Response Curves for MaxEnt Species Distribution Models
#'
#' @description
#' Computes and plots marginal response curves for each environmental variable
#' in a MaxEnt model. For each variable, a sequence of 100 values spanning its
#' observed range is evaluated while all other variables are held at their mean.
#' The result is a faceted \code{ggplot2} figure showing predicted habitat
#' suitability as a function of each predictor.
#'
#' @param maxent_mod A fitted MaxEnt model object compatible with
#'   \code{\link[dismo]{predict}} (e.g., an object of class
#'   \code{MaxEnt_model} produced by \pkg{ENMeval} 2.x).
#' @param envstack A \code{RasterStack} (\pkg{raster}) or \code{SpatRaster}
#'   (\pkg{terra}) containing the environmental variables used to fit
#'   \code{maxent_mod}. Variable names must match those used during model
#'   training.
#'
#' @return A \code{ggplot} object with one facet per environmental variable.
#'   Each panel shows predicted suitability (cloglog scale) on the y-axis
#'   and the variable's observed value range on the x-axis. The object can
#'   be further customised or saved with \code{\link[ggplot2]{ggsave}}.
#'
#' @details
#' The function implements a standard \emph{marginal response curve} approach:
#' \enumerate{
#'   \item For each variable \eqn{v}, extract the observed range from
#'     \code{envstack} and create a regular sequence of 100 equally spaced
#'     values.
#'   \item Build a prediction data frame in which \eqn{v} varies across its
#'     sequence and every other variable is fixed at its cell-level mean.
#'   \item Call \code{predict(maxent_mod, newdata, type = "cloglog")} to
#'     obtain suitability scores.
#'   \item Combine all variable-specific data frames and render a
#'     \code{facet_wrap} plot with free x-axis scales.
#' }
#' Changing the \code{type} argument inside the function body (e.g., to
#' \code{"logistic"}) will alter the suitability scale of the output curves.
#'
#' @note
#' \itemize{
#'   \item Cell values are extracted with \code{values()}, so both
#'     \pkg{raster} and \pkg{terra} objects are supported as long as the
#'     method is available for the class supplied.
#'   \item \code{NA} cells are silently removed when computing range and
#'     mean values (\code{na.rm = TRUE}).
#'   \item The function does not perform any variable importance ranking;
#'     all variables present in \code{envstack} are plotted.
#' }
#'
#' @seealso
#' \code{\link[dismo]{maxent}},
#' \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link[ggplot2]{facet_wrap}}
#'
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs theme_minimal
#'
#' @examples
#' \dontrun{
#' library(ENMeval)
#' library(terra)
#'
#' # Assume `best_mod` is a MaxEnt model selected from an ENMeval run
#' # and `env` is a SpatRaster with the training variables.
#' p <- respose_curve_function(maxent_mod = best_mod, envstack = env)
#' print(p)
#'
#' # Save to file
#' ggplot2::ggsave("response_curves.png", plot = p,
#'                 width = 12, height = 8, dpi = 300)
#' }
#'
#' @export
respose_curve_function <- function(maxent_mod, envstack) {
  
  mxmod <- maxent_mod
  envs  <- envstack
  
  vars <- names(envs)
  
  resp_list <- list()
  
  for (v in vars) {
    
    rng      <- range(values(envs[[v]]), na.rm = TRUE)
    seq_vals <- seq(rng[1], rng[2], length.out = 100)
    
    newdata           <- as.data.frame(matrix(NA, nrow = 100, ncol = length(vars)))
    colnames(newdata) <- vars
    newdata[[v]]      <- seq_vals
    
    for (vv in setdiff(vars, v)) {
      newdata[[vv]] <- mean(values(envs[[vv]]), na.rm = TRUE)
    }
    
    preds <- predict(mxmod, newdata, type = "cloglog")
    
    resp_list[[v]] <- data.frame(variable = v, x = seq_vals, y = preds)
  }
  
  resp_df <- do.call(rbind, resp_list)
  
  x <- ggplot2::ggplot(resp_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(color = "red") +
    ggplot2::facet_wrap(~ variable, scales = "free_x") +
    ggplot2::labs(x = "", y = "predict value") +
    ggplot2::theme_minimal()
  
  return(x)
}