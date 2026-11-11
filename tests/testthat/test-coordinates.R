test_that("projection reproduces the Ordnance Survey worked example", {
  # A guide to coordinate systems in Great Britain, annex C
  lat <- 52 + 39 / 60 + 27.2531 / 3600
  long <- 1 + 43 / 60 + 4.5177 / 3600
  bng <- .TransverseMercator(lat, long)
  expect_lt(abs(bng$easting - 651409.903), 1e-3)
  expect_lt(abs(bng$northing - 313177.270), 1e-3)
  ll <- .InverseTransverseMercator(651409.903, 313177.270)
  expect_lt(abs(ll$lat - lat), 1e-8)
  expect_lt(abs(ll$long - long), 1e-8)
})

test_that("LatLongToBNG() agrees with PROJ", {
  # terra::project() from EPSG:4326 to EPSG:27700, PROJ 9
  places <- data.frame(
    lat = c(54.7733, 49.9597, 60.8556, 55.2408, 52.4818),
    long = c(-1.5762, -5.2028, -0.8853, -6.5116, 1.7627),
    easting = c(427361.9, 170386.6, 460664.1, 113321.1, 655584.7),
    northing = c(542123.6, 11531.5, 1219713.1, 603351.8, 293751.7)
  )
  bng <- LatLongToBNG(places$lat, places$long)
  expect_lt(max(abs(bng$easting - places$easting)), 5)
  expect_lt(max(abs(bng$northing - places$northing)), 5)
})

test_that("BNGToLatLong() inverts LatLongToBNG()", {
  lat <- c(50.1, 54.7733, 58.5, NA)
  long <- c(-5.5, -1.5762, -3, 0)
  bng <- LatLongToBNG(lat, long)
  expect_equal(BNGToLatLong(bng$easting, bng$northing),
               data.frame(lat = lat, long = c(long[1:3], NA)),
               tolerance = 1e-6)
  expect_equal(nrow(LatLongToBNG(numeric(0), numeric(0))), 0)
})

test_that("coordinates are validated", {
  expect_error(LatLongToBNG("54", "-1"), "must be numeric")
  expect_error(LatLongToBNG(54, c(-1, -2)), "same length")
  expect_error(LatLongToBNG(91, 0), "`lat` must lie between")
  expect_error(LatLongToBNG(54, 181), "`long` must lie between")
  expect_error(BNGToLatLong(1, 1:2), "same length")
})

test_that("GridRefToBNG() reads grid references", {
  expect_equal(
    GridRefToBNG(c("NZ 274 422", "nz2742", "TG 51409 13177", "HP", NA),
                 centre = FALSE),
    data.frame(
      easting = c(427400, 427000, 651409, 400000, NA),
      northing = c(542200, 542000, 313177, 1200000, NA)
    )
  )
  expect_equal(GridRefToBNG(c("NZ 274 422", "SV")),
               data.frame(easting = c(427450, 50000),
                          northing = c(542250, 50000)))
  expect_error(GridRefToBNG(c("NZ 27", "NZ 274 42")), "\"NZ 274 42\"")
  expect_error(GridRefToBNG("IA 12 34"), "Can't read")
  expect_error(GridRefToBNG("427400"), "Can't read")
})
