#' Nearest Neighbour Distance Matching (NNDM) algorithm
#' @description
#' This function implements the \emph{NNDM} algorithm and returns the necessary
#' indices to perform a NNDM LOO CV for map validation.
#'
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param ppoints sf or sfc point object. Contains the target prediction points.
#' @param phi Numeric. Estimate of the landscape autocorrelation range in the
#' same units as the tpoints and ppoints for projected CRS, in meters for geographic CRS.
#' @param min_train Numeric between 0 and 1. Minimum proportion of training
#' data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions).
#'
#' @return An object of class \emph{nndm} consisting of a list of six elements:
#' indx_train, indx_test, and indx_exclude (indices of the observations to use as
#' training/test/excluded data in each NNDM LOO CV iteration), Gij (distances for
#' multitype G function construction between prediction and target points), Gj
#' (distances for G function construction during LOO CV), Gjstar (distances
#' for modified G function during NNDM LOO CV), phi (landscape autocorrelation range).
#' @details Details of the method can be found in the reference 'Nearest Neighbour Distance
#' Matching Leave-One-Out Cross-Validation for map validation' submitted to
#' \emph{Methods in Ecology and Evolution} (2022).
#'
#' Euclidean distances are used for projected
#' and non-defined CRS, otherwise the function uses great circle distances (units in meters).
#' @export
#' @examples
#' library("NNDM")
#' library("sf")
#'
#' # Simulate 100 random training and test points in a 100x100 square
#' set.seed(123)
#' poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
#' sample_poly <- sf::st_polygon(poly)
#' train_points <- sf::st_sample(sample_poly, 100, type = "random")
#' pred_points <- sf::st_sample(sample_poly, 100, type = "random")
#'
#' # We run NNDM. The autocorrelation range (phi) is known to be 10.
#' nndm_pred <- nndm(train_points, pred_points, 10, 0.5)
#' nndm_pred
#' plot(nndm_pred)
nndm <- function(tpoints, ppoints, phi, min_train=0){

  # If tpoints is sfc, coerce to sf.
  if(any(class(tpoints) %in% "sfc")){
    tpoints <- sf::st_sf(geom=tpoints)
  }

  # If ppoints is sfc, coerce to sf.
  if(any(class(ppoints) %in% "sfc")){
    ppoints <- sf::st_sf(geom=ppoints)
  }

  # Input data checks
  nndm_checks(tpoints, ppoints, phi, min_train)

  # Compute nearest neighbour distances between training and prediction points
  Gij <- sf::st_distance(ppoints, tpoints)
  units(Gij) <- NULL
  Gij <- apply(Gij, 1, min)

  # Compute distance matrix of training points
  tdist <- sf::st_distance(tpoints)
  units(tdist) <- NULL
  diag(tdist) <- NA
  Gj <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
  Gjstar <- Gj

  # Start algorithm
  rmin <- min(Gjstar)
  jmin <- which.min(Gjstar)[1]
  kmin <- which(tdist[jmin,]==rmin)

  while(rmin <= phi){

    # Check if removing the point improves the match. If yes, update
    if((sum(Gjstar<=rmin)-1)/length(Gjstar) >= (sum(Gij<=rmin)/length(Gij)) &
       sum(!is.na(tdist[jmin, ]))/ncol(tdist) > min_train){
      tdist[jmin, kmin] <- NA
      Gjstar <- apply(tdist, 1, function(x) min(x, na.rm=TRUE))
      rmin <- min(Gjstar[Gjstar>=rmin]) # Distances are the same for the same pair
      jmin <- which(Gjstar==rmin)[1]
      kmin <- which(tdist[jmin,]==rmin)

    }else if(sum(Gjstar>rmin)==0){
      break
    }else{ # Otherwise move on to the next distance
      rmin <- min(Gjstar[Gjstar>rmin])
      jmin <- which(Gjstar==rmin)[1]
      kmin <- which(tdist[jmin,]==rmin)
    }
  }

  # Derive indicators
  indx_train <- list()
  indx_test <- list()
  indx_exclude <- list()
  for(i in 1:nrow(tdist)){
    indx_train[[i]] <- which(!is.na(tdist[i,]))
    indx_test[[i]] <- i
    indx_exclude[[i]] <- setdiff(which(is.na(tdist[i,])), i)
  }

  # Return list of indices
  res <- list(indx_train=indx_train, indx_test=indx_test,
              indx_exclude=indx_exclude, Gij=Gij, Gj=Gj, Gjstar=Gjstar, phi=phi)
  class(res) <- c("nndm", "list")
  res
}

#' Print method for \emph{nndm} class
#'
#' @param x An object of type \emph{nndm}.
#' @param ... other arguments.
#' @export
print.nndm <- function(x, ...){
  mean_train <- round(mean(sapply(x$indx_train, length)), 2)
  min_train <- round(min(sapply(x$indx_train, length)), 2)
  cat(paste0("nndm object\n",
             "Total number of points: ", length(x$Gj), "\n",
             "Mean number of training points: ", mean_train, "\n",
             "Minimum number of training points: ", min_train, "\n"))
}

#' Plot method for \emph{nndm} class
#'
#' @param x An object of type \emph{nndm}.
#' @param ... other arguments.
#'
#' @export
plot.nndm <- function(x, ...){

  # Prepare data for plotting: Gij function
  Gij_df <- data.frame(r=x$Gij[order(x$Gij)])
  Gij_df$val <- 1:nrow(Gij_df)/nrow(Gij_df)
  Gij_df <- Gij_df[Gij_df$r <= x$phi,]
  Gij_df <- rbind(Gij_df, data.frame(r=0, val=0))
  Gij_df <- rbind(Gij_df, data.frame(r=x$phi,
                                     val=sum(x$Gij<=x$phi)/length(x$Gij)))
  Gij_df$Function <- "1_Gij(r)"

  # Prepare data for plotting: Gjstar function
  Gjstar_df <- data.frame(r=x$Gjstar[order(x$Gjstar)])
  Gjstar_df$val <- 1:nrow(Gjstar_df)/nrow(Gjstar_df)
  Gjstar_df <- Gjstar_df[Gjstar_df$r <= x$phi,]
  Gjstar_df <- rbind(Gjstar_df, data.frame(r=0, val=0))
  Gjstar_df <- rbind(Gjstar_df, data.frame(r=x$phi,
                                           val=sum(x$Gjstar<=x$phi)/length(x$Gjstar)))
  Gjstar_df$Function <- "2_Gjstar(r)"

  # Prepare data for plotting: G function
  Gj_df <- data.frame(r=x$Gj[order(x$Gj)])
  Gj_df$val <- 1:nrow(Gj_df)/nrow(Gj_df)
  Gj_df <- Gj_df[Gj_df$r <= x$phi,]
  Gj_df <- rbind(Gj_df, data.frame(r=0, val=0))
  Gj_df <- rbind(Gj_df, data.frame(r=x$phi,
                                   val=sum(x$Gj<=x$phi)/length(x$Gj)))
  Gj_df$Function <- "3_Gj(r)"

  # Merge data for plotting
  Gplot <- rbind(Gij_df, Gjstar_df, Gj_df)

  # Plot
  ggplot2::ggplot(Gplot) +
    ggplot2::geom_step(ggplot2::aes_string(x="r", y="val", colour="Function", size="Function"),
                       alpha = 0.8) +
    ggplot2::scale_size_manual(values=c(1.1, 1.1, 0.5),
                               labels=c(expression(hat(G)[ij](r)),
                                        expression(hat(G)[j]^"*"*"(r,"*bold(L)*")"),
                                        expression(hat(G)[j](r)))) +
    ggplot2::scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9"),
                                 labels=c(expression(hat(G)[ij](r)),
                                          expression(hat(G)[j]^"*"*"(r,"*bold(L)*")"),
                                          expression(hat(G)[j](r)))) +
    ggplot2::ylab(expression(paste(hat(G)[ij](r), ", ",
                                   hat(G)[j]^"*"*"(r,"*bold(L)*")", ", ",
                                   hat(G)[j](r)))) +
    ggplot2::labs(colour="", size="") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.text.align=0,
                   legend.text=ggplot2::element_text(size=12))
}

#' Input data checks for NNDM
#'
#' @param tpoints sf or sfc point object. Contains the training points samples.
#' @param ppoints sf or sfc point object. Contains the prediction target points.
#' @param phi Numeric. Estimate of the landscape autocorrelation range (in the
#' same units as the tpoints and ppoints).
#' @param min_train Numeric between 0 and 1. Minimum proportion of training
#' data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions).
nndm_checks <- function(tpoints, ppoints, phi, min_train){

  # Check for valid range of phi
  if(phi < 0){
    stop("phi must be positive.")
  }

  # min_train must be a single positive integer
  if(length(min_train)!=1 | min_train<0 | min_train>1 | !is.numeric(min_train)){
    stop("min_train must be a numeric between 0 and 1.")
  }

  # Check class and geometry type of tpoints
  if(!any(c("sfc", "sf") %in% class(tpoints))){
    stop("tpoints must be a sf/sfc object.")
  }else if(!any(class(sf::st_geometry(tpoints)) %in% c("sfc_POINT"))){
    stop("tpoints must be a sf/sfc point object.")
  }

  # Check class and geometry type of ppoints
  if(!any(c("sfc", "sf") %in% class(ppoints))){
    stop("ppoints must be a sf/sfc object.")
  }else if(!any(class(sf::st_geometry(ppoints)) %in% c("sfc_POINT"))){
    stop("ppoints must be a sf/sfc point object.")
  }

  # Check same CRS of tpoints and ppoints
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
    stop("tpoints and ppoints must have the same CRS.")
  }
}
