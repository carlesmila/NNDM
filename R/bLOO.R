#' Creates indices for bLOO CV
#' @description
#' Function that computes the indices necessary for bLOO CV given a set of training
#' points and a buffer radius.
#' @param tpoints sf/sfc point object. Contains the training points samples.
#' @param radius Numeric. Positive radius of the exclusion buffer  (in the same units as the tpoints).
#' @param min_train Numeric, between 0 and 1. Minimum proportion of training
#' data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions). If
#' positive and not met for a given iteration, the radius will be decreased by
#' 1/10 until the condition is met.
#' @return An object of class \emph{bLOO} consisting of a list of four elements:
#' indx_train, indx_test, and indx_exclude (indices of the observations to use as
#' training/test/excluded data in each bLOO CV fold), and radii (useful if min_train>0,
#' as the radius indicated by the user might have been reduced for a given iteration).
#' @export
#' @details
#' Euclidean distances are used for projected and non-defined CRS, otherwise the
#' function uses great circle distances.
#' @examples
#' library("NNDM")
#' library("sf")
#'
#' # Simulate 100 random points in a 100x100 square
#' set.seed(1234)
#' poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' sample_poly <- sf::st_polygon(poly)
#' sim_points <- sf::st_sample(sample_poly, 100, type = "random")
#'
#' # Compute bLOO with a radius of 40 and no training data size restrictions.
#' bLOO40 <- bLOO(sim_points, 40)
#' bLOO40
#'
#' # Visualize results for the first iteration
#' plot(bLOO40, 1)
#'
#' # Compute bLOO with a radius of 40 and a minimum training data size of 50%.
#' bLOO40_mintrain <- bLOO(sim_points, 40, 0.5)
#' bLOO40_mintrain
bLOO <- function(tpoints, radius, min_train=0){

  # Input data checks
  bLOO_checks(tpoints=tpoints, radius=radius, min_train=min_train)

  # If tpoints is sfc, coerce to sf.
  if(any(class(tpoints) %in% "sfc")){
    tpoints <- sf::st_sf(geom=tpoints)
  }

  # Create empty lists to store indices
  indx_train <- list()
  indx_test <- list()
  indx_exclude <- list()
  radii <- list()

  # Create indicators for referencing
  tpoints$ind <- 1:nrow(tpoints)

  # Compute distance matrix
  dists <- sf::st_distance(tpoints)

  # Start looping through all the tpoints
  for(i in 1:nrow(tpoints)){

    # Create subsets of data and compute train and excluded data
    disti <- as.numeric(dists[,i])
    exc_tpoints <- tpoints[disti<=radius & disti!=0,]
    tr_tpoints <- tpoints[disti>radius,]

    # Reduce the radius to meet minimum training size requirements?
    j <- 1
    if(min_train>0 & nrow(tr_tpoints)/nrow(tpoints)<min_train){
      while(nrow(tr_tpoints)/nrow(tpoints)<min_train){
        j <- j-0.1
        exc_tpoints <- tpoints[disti<=radius*j & disti!=0,]
        tr_tpoints <- tpoints[disti>radius*j,]
      }
    }

    # Store indices and radius
    indx_train[[i]] <- tr_tpoints$ind
    indx_test[[i]] <- i
    indx_exclude[[i]] <- exc_tpoints$ind
    radii[[i]] <- radius*j
  }

  # Return list of indices
  res <- list(indx_train=indx_train, indx_test=indx_test,
              indx_exclude=indx_exclude, radii=radii)
  class(res) <- c("bLOO", "list")
  attr(res, "data") <- tpoints
  res
}


#' Print method for \emph{bLOO} class
#'
#' @param x An object of type \emph{bLOO}.
#' @param ... other arguments.
#' @export
print.bLOO <- function(x, ...){
  cat(paste0("bLOO object\n",
             "Total number of points: ", nrow(attr(x, "data")), "\n",
             "Mean number of training points: ",
             round(mean(sapply(x$indx_train, length)), 2), "\n",
             "Minimum number of training points: ",
             round(min(sapply(x$indx_train, length)), 2), "\n",
             "Mean buffer radius: ",
             round(mean(unlist(x$radii, length)),2), "\n"))
}


#' Plot method for \emph{bLOO} class
#'
#' @param x An object of type \emph{bLOO}.
#' @param it bLOO iteration that needs to be plotted.
#' @param ... other arguments.
#'
#' @export
plot.bLOO <- function(x, it, ...){

  # Ensure it is lower than the number of training points
  if(it > length(x$indx_train)){
    stop(paste0("The bLOO only has ", length(x$indx_train),
                " iterations"))
  }

  # Prepare data for plotting
  plotdata <- attr(x, "data")
  plotdata$role <-  NA
  plotdata$role[x$indx_test[[it]]] <- "test"
  plotdata$role[x$indx_train[[it]]] <- "train"
  plotdata$role[x$indx_exclude[[it]]] <- "exclude"

  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_sf(data=plotdata, ggplot2::aes_string(colour="role")) +
    ggplot2::scale_colour_viridis_d(end=0.9) +
    ggplot2::ggtitle(paste0("bLOO CV iteration ", it)) +
    ggplot2::theme_bw()
}


#' Input data checks for bLOO
#'
#' @param tpoints sf point object. Contains the training points samples.
#' @param radius Numeric. Positive radius of the exclusion buffer (in the same units as the tpoints).
#' @param min_train Numeric, between 0 and 1. Proportion of minimum training
#' data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions). If
#' positive and not met for a given iteration, the radius will be decreased by
#' 1/10 until the condition is met.
bLOO_checks <- function(tpoints, radius, min_train){
  # Check valid range of radius and min_train
  if(radius < 0){
    stop("Radius must be positive.")
  }
  if(min_train>1|min_train<0){
    stop("min_train must be between 0 and 1.")
  }
  # Check that input is a sf or sfc point object.
  if(!any(class(tpoints) %in% c("sf", "sfc"))){
    stop("Input tpoints must be of class sf or sfc.")
  }
  if(!any(class(sf::st_geometry(tpoints)) %in% c("sfc_POINT"))){
    stop("Wrong type of geometry of input data.")
  }
}
