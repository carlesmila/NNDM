% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nndm.R
\name{nndm}
\alias{nndm}
\title{Nearest Neighbour Distance Matching (NNDM) algorithm}
\usage{
nndm(tpoints, ppoints, phi, min_train = 0)
}
\arguments{
\item{tpoints}{sf or sfc point object. Contains the training points samples.}

\item{ppoints}{sf or sfc point object. Contains the target prediction points.}

\item{phi}{Numeric. Estimate of the landscape autocorrelation range in the
same units as the tpoints and ppoints for projected CRS, in meters for geographic CRS.}

\item{min_train}{Numeric between 0 and 1. Minimum proportion of training
data that must be used in each CV fold. Defaults to 0 (i.e. no restrictions).}
}
\value{
An object of class \emph{nndm} consisting of a list of six elements:
indx_train, indx_test, and indx_exclude (indices of the observations to use as
training/test/excluded data in each NNDM LOO CV iteration), Gij (distances for
multitype G function construction between prediction and target points), Gj
(distances for G function construction during LOO CV), Gjstar (distances
for modified G function during NNDM LOO CV), phi (landscape autocorrelation range).
}
\description{
This function implements the \emph{NNDM} algorithm and returns the necessary
indices to perform a NNDM LOO CV for map validation.
}
\details{
Details of the method can be found in the reference 'Nearest Neighbour Distance
Matching Leave-One-Out Cross-Validation for map validation' submitted to
\emph{Methods in Ecology and Evolution} (2022).

Euclidean distances are used for projected
and non-defined CRS, otherwise the function uses great circle distances (units in meters).
}
\examples{
library("NNDM")
library("sf")

# Simulate 100 random training and test points in a 100x100 square
set.seed(123)
poly <- list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2, byrow=TRUE))
sample_poly <- sf::st_polygon(poly)
train_points <- sf::st_sample(sample_poly, 100, type = "random")
pred_points <- sf::st_sample(sample_poly, 100, type = "random")

# We run NNDM. The autocorrelation range (phi) is known to be 10.
nndm_pred <- nndm(train_points, pred_points, 10, 0.5)
nndm_pred
plot(nndm_pred)
}
