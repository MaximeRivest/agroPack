#' Run linear mixed effect models on the soil properties as a function of windbreaks
#'
#' This function plot the soil property given for example:soil water content,
#' water holdin capacity, organic C, total N, C:N ratio, and soil pH and regress it to the distance
#' from the windbreak while having fields as a random effect.
#'
#' @usage stats_for_soil_properties(y = main_data$Whc)
#' @param y a vector of data from the main_data dataset.
#' @details No details for now.
#' @return list of relevant stats.
#' @author Maxime Rivest
#' @examples
#' stats <- stats_for_soil_properties(y = main_data$Whc)
#' @seealso \code{\link{agroPack}}.
#' @keywords manip
#' @export
#' @import ggplot2
plot_for_soil_properties <- function(y, my_ylab = ""){
  df <- y$df
  df$Treatment_meter <- as.numeric(df$Treatment)
  df$Treatment_meter[df$Treatment_meter==1] <- 0
  df$Treatment_meter[df$Treatment_meter==2] <- 8
  df$Treatment_meter[df$Treatment_meter==3] <- 50
  df$Treatment_meter <- as.factor(df$Treatment_meter)

  # Annotate for significance
  p12 <- y$summary_model123$coefficients[2,5]
  p13 <- y$summary_model123$coefficients[3,5]
  p23 <- y$summary_model231$coefficients[2,5]

  if(p12 < 0.05 & p13 < 0.05 & p23 < 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","b","c"))
  } else if(p12 < 0.05 & p13 < 0.05 & p23 >= 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","b","b"))
  } else if(p12 < 0.05 & p13 >= 0.05 & p23 < 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","b","a"))
  } else if(p12 >= 0.05 & p13 < 0.05 & p23 < 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","a","b"))
  } else if(p12 < 0.05 & p13 >= 0.05 & p23 >= 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","b","ab"))
  } else if(p12 >= 0.05 & p13 >= 0.05 & p23 < 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("ab","a","b"))
  } else if(p12 >= 0.05 & p13 < 0.05 & p23 >= 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y),
                          texti = c("a","ab","b"))
  } else if(p12 >= 0.05 & p13 >= 0.05 & p23 >= 0.05) {
    anno_df <- data.frame(x = c(1,2,3),
                          y = max(df$y)+0.1*max(df$y) ,
                          texti = c("a","a","a"))
  }

  my_plot <- ggplot2::ggplot(df,
                             ggplot2::aes(x = Treatment_meter,
                                          y = y
                                          #,color = as.factor(Treatment_meter)
                             ))+
    ggplot2::geom_boxplot(alpha = 1)+
    ggplot2::xlab("Distance from windbreak (m)") +
    ggplot2::ylab(my_ylab)+
    theme_Publication() +
    ggplot2::stat_summary(fun.y=mean,
                          colour="black",
                          geom="point",
                          shape=1,
                          size=2,
                          show.legend = FALSE) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::geom_text(data = anno_df, mapping = ggplot2::aes(x, y, label = texti))

  return(my_plot)
}
