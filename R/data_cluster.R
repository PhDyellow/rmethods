#' Helper function to allow both warnings and values to be captured
#'
#' @param expr a code block to evaluate
#'
#' @return list of $value and $warnings, $value contains the result of evaluating the code block regardless of warnings,
#' $warnings contains a list of warnings generated during evaluation



#https://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}


#' Find optimal clusters by according to ICL
#'
#' @param data A matrix or data.frame to cluster
#' @param mclust_model Mclust model string such as "EII". see ?mclustModelNames
#' @param n_groups Vector of group counts to test
#' @param n_iter Number of times to refit each group count in n_groups
#' @param min_site_ratio Only accept fits with all clusters owning at least \code{nrow(data) * min_site_ratio} sites
#'
#' @return list of scores and best mclust object
#'
#' @importFrom foreach foreach %:%
#' @import mclust

optimal_icl <- function(data, mclust_model, n_groups, n_iter = 3, min_site_ratio){


#to reduce the memory footprint, only the best model fitted so far is kept.
#each parallel worker should keep local best models, global best will be evaluated later
max_icl <- NULL
max_model <- NULL

if(!foreach::getDoParRegistered()){
  foreach::registerDoSEQ()
}

map_wrap <- function(...){mapply(c, SIMPLIFY = FALSE, ...)}

scores <- foreach::foreach(g = n_groups, .combine = map_wrap) %:%
  foreach::foreach(i = 1:n_iter, .combine = map_wrap, .packages = c("mclust", "rmethods")) %dopar% {

    icl_test <- rmethods::withWarnings(mclust::mclustICL(data, G = g, modelNames = mclust_model, warn = TRUE))

    if(!is.null(icl_test$warnings)){
      #Warning in tryCatch loop, usually empty cluster, sometimes fitting singularity
      #I will keep "no assignment" errors, but not others
      if(length(icl_test$warnings) > 1 | grepl(pattern = "^no assignment", x = icl_test$warnings[[1]]$message)){
        message(paste0("Discarded ", g, ", because:", icl_test$warnings))
        return(list(groups = g, icl = max(icl_test$value), type = "empty", fit = list(NA)))
      }

    }


    bic_test <- icl_test$value
    class(bic_test) <- "mclustBIC" #hack to make ICL as functional as BIC in mclust
    mclust_fit <- rmethods::withWarnings(mclust::Mclust(data, x = bic_test))

    if(!is.null(mclust_fit$warnings)){
      #Warning in tryCatch loop, usually empty cluster, sometimes fitting singularity
      if(length(icl_test$warnings) > 1 | grepl(pattern = "^no assignment", x = icl_test$warnings[[1]]$message)){
        message(paste0("Discarded ", g, ", because:", mclust_fit$warnings))
        return(list(groups = g, icl = max(icl_test$value), type = "empty2", fit = list(NA)))
      }
    }
    mclust_fit <- mclust_fit$value

    if (min(table(mclust_fit$classification)) < nrow(data) * min_site_ratio){
      #A site has less than min_site_ratio sites
      message(paste0("Rejected ", g, ". A cluster had ", min(table(mclust_fit$classification)),  ", less than ", nrow(data) * min_site_ratio, " sites.\n"))
      return(list(groups = g, icl = max(icl_test$value), type = "empty", fit = list(NA)))
    } else if(is.null(max_icl)){
      #The first good model has been found,
      #Model returned no errors and has at least min_site_ratio sites in all clusters
      max_icl <- icl_test$value
      max_model <- mclust_fit
      message(paste0("Accepted ", g, ". First valid fitting.\n"))
      return(list(groups = g, icl = max(icl_test$value), type = "peak", fit = list(max_model)))
    } else if(max(icl_test$value) > max(max_icl)){
      message(paste0("Accepted ", g, ". ", max(icl_test$value), " > ", max(max_icl), "\n"))
      max_icl <- icl_test$value
      max_model <- mclust_fit
      return(list(groups = g, icl = max(icl_test$value), type = "peak", fit = list(max_model)))
    } else {
      message(paste0("Did not accept ", g, ". Below maximum ICL.\n"))
      return(list(groups = g, icl = max(icl_test$value), type = "below", fit = list(NA)))
    }
    message("Hit default")
    return(list(groups = g, icl = max(icl_test$value), type = "default", fit = list(NA)))
  }


valid <- scores$type == "peak"
bestICL <- scores$fit[valid][[which.max(scores$icl[valid])]]
return(list(scores = data.frame(scores[c("groups", "icl", "type")]), best_fit = bestICL))

}
