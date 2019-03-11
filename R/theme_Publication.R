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
theme_Publication <- function(base_size=12, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = margin(unit(0, "cm")),
            legend.title = element_blank(),
            legend.background = element_rect(colour = NA),
            legend.box.background = element_blank(),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text()
    ))
}
