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
    message(paste0("Fitting, G = ", g, ", iteration ", i))

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

#' Cluster by connectivity, greedy
#'
#' Generates a clustering from a binary connectivity matrix.
#'
#' The first cluster center is the site with the most connected neighbours.
#' The center and all center neighbours are the first cluster.
#' The second cluster center is the site with the most connected neighbours excluding any sites in the first cluster.
#' Each new cluster is formed from sites not yet connected to any center.
#'
#' @param conn_mat binary site by site connectivity matrix, TRUE/1 if sites i and j are connected
#' @param stop_all Should sites with no neighbours be returned as cluster centers? If TRUE, yes, if FALSE, no.
#' @param break_ties How should ties be broken when multiple sites tie for having the most connected neighbours. One of "random", "first"
#' @param min_sites If the site with the most neighbours has less than min_sites neighbours, do not return any more cluster centers
#' @param max_primary Maximum number of clusters to allow.
#'
#' @return vector of site row indicies that are cluster centers.
#'
#' @examples
#' Membership in each cluster is then given by conn_mat[, centers]
#'
#' binary_membership <- (rowSums((membership*1) * matrix(2^(seq(ncol(membership), 1)), nrow = nrow(membership), ncol = ncol(membership), byrow = TRUE)))
#' binary_membership_log <- log2(binary_membership)
conn_clust_rec <- function(conn_mat, stop_all = TRUE, break_ties = "random", min_sites = 100, max_primary = Inf){
  #Function to assign sites into clusters
  # Given a connectivity matrix, where 1/TRUE indictates connected and not significantly different
  #Returns the row indexes of cluster centers, sorted from largest to smallest clusters
  # the cluster membership can be found by conn_mat[,ind]
  #Diagonals must be all 1/TRUE

  if(any(diag(conn_mat) != 1)){
    stop("Connectivity matrix must have 1/TRUE on the diagonals")
  }

  #Function is recursive.
  #Given a conn_mat, finds the next largest centre, removes all sites connected to the centre,
  #and passes down the subset conn_mat to the next level.

  #Stopping condition: All sites have connectivity only to themselves.
  conn_total <- rowSums(conn_mat)
  conn_total_max <- max(conn_total)

  if(conn_total_max ==1 ){
    #All sites are connected only to themselves.
    if(stop_all){
      #Return all rows as centers
      return(1:nrow(conn_mat))
    } else {
      #Return NO rows as centers
      #Dont return any rows. Some sites may end up having NO cluster membership, check later
      return(numeric(0))
    }
  }


  #Not at stopping condition, pick a center
  #all maxima
  conn_max_ind <- which(conn_total %in% max(conn_total))

  new_center <- switch(EXPR = break_ties,
                       random = {conn_max_ind[sample.int(length(conn_max_ind), 1)]},
                       first = {conn_max_ind[1]},
                       default = {stop("break_ties not specified properly")})
  #the connectivity matrix is a logical matrix of sites connected to sites
  #If I take the row of new_centre, then keep only sites that are FALSE, then the
  #remaining sub_conn_mat is all sites not connected to the new_centre
  connected_sites <- conn_mat[new_center, ]
  sub_conn_mat <- conn_mat[!connected_sites , !connected_sites]

  if(dim(sub_conn_mat)[1] == 0 & dim(sub_conn_mat)[2] == 0) {
    #All sites are connected to new_centre
    return(new_center)
  }

  if(dim(sub_conn_mat)[1] < min_sites & dim(sub_conn_mat)[2] < min_sites ) {
    #Too few sites left. Don't go deeper
    return(new_center)
  }

  if(max_primary == 1) {
    #Reached max depth
    return(new_center)
  }

  sub_center_ind <- conn_clust_rec(sub_conn_mat, stop_all = stop_all, break_ties = break_ties, min_sites = min_sites, max_primary = max_primary -  1)

  #sub_ind contains all the centers from the sub_conn_mat.
  #translate into current conn_mat row indexes, then return current new_center at front of
  #all row_indexes

  sub_site_ind <- 1:sum(!connected_sites)

  translate_vec <- rep(0, length.out = length(connected_sites))
  translate_vec[!connected_sites] <- sub_site_ind

  sub_center_current <- foreach(center = sub_center_ind, .inorder = TRUE, .combine = c) %do% {
    which(translate_vec == center)
  }

  return(c(new_center, sub_center_current))
}
