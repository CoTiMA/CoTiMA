#' ctmaCompFit
#'
#' @param model1 ""
#' @param model2 ""
#'
#' @return
#' @export
#'
#'
ctmaCompFit <- function(model1=NULL, model2=NULL) {
  n.model1 <- length(model1$summary$model); n.model1
  n.model2 <- length(model2$summary$model); n.model2
  names.model1 <- model1$summary$model; names.model1
  names.model2 <- model2$summary$model; names.model2
  allResults <- list()
  counter <- 0
  for (i in 1:n.model1) {
    for (j in 1:n.model2) {
      counter <- counter + 1; counter
      allResults[[counter]] <- list()
      if (counter > 1) {
        allResults[[counter]]$" " <- "############################    MODEL COMPARISON   ############################"
      } else {
          allResults[[counter]]$" " <- "##########################   NEXT MODEL COMPARISON   ##########################"
        }
      allResults[[counter]]$"Model: " <- paste0(names.model1[i], " COMPARED WITH   Model: ", names.model2[j]); allResults[[counter]]$comparison
      allResults[[counter]]$"Diff_Minus2LL: " <- abs((abs(model1$summary$minus2ll[[i]]) - abs(model2$summary$minus2ll[[j]]))); allResults[[counter]]$Diff_Minus2LL
      allResults[[counter]]$"Diff_df: " <- abs((model1$summary$n.parameters[[i]] - model2$summary$n.parameters[[j]])); allResults[[counter]]$Diff_df
      if (allResults[[counter]]$"Diff_df: " > 0 ) {
        allResults[[counter]]$"prob: " <- 1-stats::pchisq(allResults[[counter]]$Diff_Minus2LL, allResults[[counter]]$Diff_df); allResults[[counter]]$prob
      } else {
        allResults[[counter]]$"prob: " <- "Diff_df = 0, thus probability cannot be computed"
      }
      allResults[[counter]]$"Message1: " <- "A prob value < .05 indicates a significant difference."; allResults[[counter]]$Message
      allResults[[counter]]$"Message2: " <- "This is valid only if models are nested, which is NOT checked"; allResults[[counter]]$Message
    }
  }
  return(list(paste(names(unlist(allResults)), unlist(allResults))))
}

