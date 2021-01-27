#' ctmaOptimzeInit
#'
#' @param oldStudyList oldStudyList
#' @param activeDirectory activeDirectory
#' @param problemStudy problemStudy
#' @param reFits reFits
#' @param n.latent n.latent
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#'
#' @export ctmaOptimzeInit
#'
ctmaOptimzeInit <- function(oldStudyList=NULL,
                            activeDirectory=NULL,
                            problemStudy=NULL,
                            reFits=NULL,
                            n.latent=NULL)
  {
  doParallel::registerDoParallel(detectCores()-1)
  '%dopar%' <- foreach::'%dopar%'

    if (is.null(oldStudyList) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop("arguments are missing")
    newStudyList <- list()
    for (j in 1:(length(names(oldStudyList))-1)) {
      tmp1 <- lapply(oldStudyList[[names(oldStudyList)[j]]], function(x) x)
      newStudyList[[j]] <- list(tmp1[[problemStudy]])
    }
    newStudyList[[length(names(oldStudyList))]] <- 1
    names(newStudyList) <- names(oldStudyList)
    allfits <- foreach::foreach(i=1:reFits) %dopar% {
      fits <- ctmaInit(newStudyList, coresToUse = 1, n.latent=n.latent,
                       activeDirectory = activeDirectory)
      return(fits)
    }
    all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
    bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
    bestFit <- allfits[[bestFit]]
    invisible(list(bestFit=bestFit, all_minus2ll=all_minus2ll))
  }


