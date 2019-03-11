#' Present the results of a mixed effect models on the soil properties as a function of windbreaks
#'
#' Generate a clean data frame.
#'
#' @usage table_of_multiple_comparaison(my_list, varname = "Name to display")
#' @param my_list a list out of the function stats_for_soil_properties
#' @param varname a character vector of length one holding the name of the variable that should
#'  be displayed in the table
#' @details No details for now.
#' @return data.frame of coefficients for multiple comparaisons
#' @author Maxime Rivest
#' @examples
#' stats <- stats_for_soil_properties(y = main_data$Whc)
#' table_of_multiple_comparaison(my_list = stats, varname = "Water holding capacity")
#' @seealso \code{\link{agroPack}}.
#' @keywords manip
#' @export
table_of_multiple_comparaison <- function(my_list, varname = ""){
  coef123 <- as.data.frame(my_list$summary_model123$coefficients)
  #coef123$Estimate <- coef123$Estimate + coef123$Estimate[1]
  coef123 <- coef123[2:3,]
  coef123$Comparaison <- c("0m to 8m", "0m to 50m")
  coef231 <- as.data.frame(my_list$summary_model231$coefficients)
  #coef231$Estimate <- coef231$Estimate + coef231$Estimate[1]
  coef231 <- coef231[2,]
  coef231$Comparaison <- c("8m to 50m")
  coefs <- rbind(coef123,coef231)
  row.names(coefs) <- 1:3
  coefs <- cbind(data.frame("soil property" = rep(varname, 3)),coefs$Comparaison, coefs[,1:5])
  names(coefs) <- c("soil property", "pair-wise comparaison", "difference between means", "Standard error", "degrees of freedom", "t-value", "p-value")
  return(coefs)
}
