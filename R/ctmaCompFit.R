#' ctmaCompFit
#'
#' @description Performs log-liklihood ratio tests to compare the fit of 2 models (CoTiMAFit objects created with \code{\link{ctmaFit}}
#' or \code{\link{ctmaEqual}}), i.e., the difference between the two -2 times LLs between the first model and the more constrained second
#' model. The nested structure of the two models is assumed to be given and not checked.
#'
#' @param model1 Model 1
#' @param model2 Model 2
#'
#' @importFrom stats pchisq
#'
#' @examples
#' minus2llDiffTest <- ctmaCompFit(CoTiMAFullInv23Fit_6,
#'                                 CoTiMAFullInvEq23Fit_6)
#' summary(minus2llDiffTest)
#'
#' @export ctmaCompFit
#'
#' @return Returns the the difference between the two -2 times LLs (Diff_Minus2LL), the associated difference in degrees of
#' freedom (Diff_df (= Diff_n.params)), and the probability (prob).
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
      allResults[[counter]]$"Diff_df (= Diff_n.params): " <- abs((model1$summary$n.parameters[[i]] - model2$summary$n.parameters[[j]])); allResults[[counter]]$Diff_df
      if (allResults[[counter]]$"Diff_df (= Diff_n.params): " > 0 ) {
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

