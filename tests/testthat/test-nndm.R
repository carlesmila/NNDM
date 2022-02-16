# Compute data for testing
set.seed(1234)
poly <- sf::st_polygon(list(matrix(c(0,0,0,50,50,50,50,0,0,0), ncol=2,
                                   byrow=TRUE)))
tpoints_sfc <- sf::st_sample(poly, 50, type = "random")
tpoints_sf <- sf::st_sf(geom = tpoints_sfc)
ppoints_sfc <- sf::st_sample(poly, 50, type = "regular")
ppoints_sf <- sf::st_sf(geom = ppoints_sfc)

# Mercator: 3857
# WGS84: 4326

# NNDM testing
test_that("Valid range of phi", {
  expect_error(nndm(tpoints_sf, ppoints_sf, -1, 0.5),
               "phi must be positive.")
})

test_that("NNDM detects wrong data and geometry types", {
  # tpoints
  expect_error(nndm(1, ppoints_sf, 10, 0.5),
               "tpoints must be a sf/sfc object.")
  expect_error(nndm(poly, ppoints_sf, 10, 0.5),
               "tpoints must be a sf/sfc object.")
  expect_error(nndm(sf::st_sfc(poly), ppoints_sf, 10, 0.5),
               "tpoints must be a sf/sfc point object.")
  # ppoints
  expect_error(nndm(tpoints_sf, 1, 10, 0.5),
               "ppoints must be a sf/sfc object.")
  expect_error(nndm(tpoints_sf, poly, 10, 0.5),
               "ppoints must be a sf/sfc object.")
  expect_error(nndm(tpoints_sf, sf::st_sfc(poly), 10, 0.5),
               "ppoints must be a sf/sfc point object.")
})

test_that("NNDM detects different CRS in inputs", {

  tpoints_sf_4326 <- sf::st_set_crs(tpoints_sf, 4326)
  tpoints_sf_3857 <- sf::st_set_crs(tpoints_sf, 3857)
  ppoints_sf_4326 <- sf::st_set_crs(ppoints_sf, 4326)
  ppoints_sf_3857 <- sf::st_set_crs(ppoints_sf, 3857)

  # tests
  expect_error(nndm(tpoints_sf_3857, ppoints_sf, 10, 0.5),
               "tpoints and ppoints must have the same CRS.")
  expect_error(nndm(tpoints_sf_3857, ppoints_sf_4326, 10, 0.5),
               "tpoints and ppoints must have the same CRS.")
})

test_that("NNDM yields the expected results for all data types in the docs", {

  # tpoints, ppoints, and sarea as sf
  expect_equal(as.numeric(round(
    nndm(tpoints_sf, ppoints_sf, 10, 0.5)$Gjstar[1], 4)), 3.7266)
  # tpoints as sfc
  expect_equal(as.numeric(round(
    nndm(tpoints_sfc, ppoints_sf, 10, 0.5)$Gjstar[2], 4)), 2.7805)
  # ppoints as sfc
  expect_equal(as.numeric(round(
    nndm(tpoints_sf, ppoints_sfc, 10, 0.5)$Gjstar[3], 4)), 8.0966)
  # sarea as sfc
  expect_equal(as.numeric(round(
    nndm(tpoints_sf, ppoints_sf, 10, 0.5)$Gjstar[4], 4)), 6.0793)
  # sarea as sfg
  expect_equal(as.numeric(round(
    nndm(tpoints_sf, ppoints_sf, 10, 1)$Gjstar[5], 4)), 4.9418)
  # vary ratio
  expect_equal(as.numeric(round(
    nndm(tpoints_sf, ppoints_sf, 10, 0.5)$Gjstar[5], 4)), 4.9418)
})

test_that("NNDM yields the expected results for all CRS", {

  # Projected
  tpoints_3857 <- sf::st_set_crs(tpoints_sf, 3857)
  ppoints_3857 <- sf::st_set_crs(ppoints_sf, 3857)
  expect_equal(as.numeric(round(
    nndm(tpoints_3857, ppoints_3857, 10, 0.5)$Gjstar[1], 4)), 3.7266)

  # Geographic
  tpoints_sf_4326 <- sf::st_set_crs(tpoints_sf, 4326)
  ppoints_sf_4326 <- sf::st_set_crs(ppoints_sf, 4326)
  expect_equal(as.numeric(round(
    nndm(tpoints_sf_4326, ppoints_sf_4326, 10, 0.5)$Gjstar[5], 4)), 546757.98)
})
