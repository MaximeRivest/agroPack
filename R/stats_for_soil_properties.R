#' Run linear mixed effect models on the soil properties as a function of windbreaks
#'
#' This function statistically tests the soil property given for example:soil water content,
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
#' @import lmerTest
stats_for_soil_properties <- function(y){
  df <- data.frame(Field_block = main_data$Field_block,
                   Treatment = main_data$Treatment,
                   y = y)
  df$Field_block <- as.factor(df$Field_block)
  df$Treatment <- as.factor(df$Treatment)
  df <- na.omit(df)
  mean_y <- mean(y, na.rm = T)
  sd_y <- sd(y, na.rm = T)

  my_model123 <- lmerTest::lmer(y ~ Treatment + (1|Field_block),data = df)

  reorderdf <- df
  reorderdf$Treatment <- factor(reorderdf$Treatment,levels(reorderdf$Treatment)[c(2,3,1)])
  my_model231 <- lmerTest::lmer(y ~ Treatment + (1|Field_block),data = reorderdf)

  list_of_res <- list("df" = df,
                      "mean_y" = mean_y,
                      "sd_y" = sd_y,
                      "summary_model123" = summary(my_model123),
                      "summary_model231" = summary(my_model231))
  return(list_of_res)
}
