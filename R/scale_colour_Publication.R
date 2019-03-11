#' Ggplot theme for publication
#'
#' A theme
#'
#' @usage theme_Publication()
#' @param y a vector of data from the main_data dataset.
#' @details No details for now.
#' @return list of relevant stats.
#' @author Maxime Rivest
#' @examples
#' stats <- stats_for_soil_properties(y = main_data$Whc)
#' @seealso \code{\link{agroPack}}.
#' @keywords manip
#' @export
#' @import ggplot2 grid ggthemes
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
