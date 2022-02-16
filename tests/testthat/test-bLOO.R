# Compute data for testing
set.seed(1234)
poly <- sf::st_polygon(list(matrix(c(0,0,0,100,100,100,100,0,0,0), ncol=2,
                                    byrow=TRUE)))
spoints <- sf::st_sample(poly, 100, type = "random")

# Start tests
test_that("function detects wrong radii", {
  expect_error(bLOO(spoints, -1), "Radius must be positive.")
})

test_that("function detects wrong min_train", {
  expect_error(bLOO(spoints, 40, 3), "min_train must be between 0 and 1.")
})

test_that("function detects wrong data types", {
  expect_error(bLOO(1:10, 40), "Input tpoints must be of class sf or sfc.")
})

test_that("function detects wrong geometry types", {
  expect_error(bLOO(sf::st_sfc(poly), 40),
               "Wrong type of geometry of input data.")
})

test_that("function can deal with sfc", {
  restest <- bLOO(spoints, 40)
  expect_equal(sapply(restest$indx_train, length)[1], 80)
  expect_equal(sapply(restest$indx_exclude, length)[1], 19)
  expect_equal(unique(sapply(restest$indx_test, length)), 1)
  expect_equal(unique(unlist(restest$radii)), 40)
})

test_that("function can deal with sf", {
  restest <- bLOO(sf::st_sf(geom=spoints), 40)
  expect_equal(sapply(restest$indx_train, length)[2], 55)
  expect_equal(sapply(restest$indx_exclude, length)[2], 44)
  expect_equal(unique(sapply(restest$indx_test, length)), 1)
  expect_equal(unique(unlist(restest$radii)), 40)
})

test_that("results with min_train!=0", {
  restest <- bLOO(spoints, 40, min_train=0.6)
  expect_equal(sapply(restest$indx_train, length)[38], 65)
  expect_equal(sapply(restest$indx_exclude, length)[38], 34)
  expect_equal(unique(sapply(restest$indx_test, length)), 1)
  expect_equal(unique(unlist(restest$radii)), c(40,36,32))
})

test_that("plot detects wrong iteration in input", {
  expect_error(plot(bLOO(spoints, 40), 130), "The bLOO only has ")
})
