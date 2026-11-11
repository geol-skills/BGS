test_that("lookups agree whichever way a point is given", {
  # Durham Cathedral
  byLatLong <- Bedrock(54.7733, -1.5762)
  expect_equal(byLatLong$lex, "PSMCM")
  expect_equal(byLatLong$maxPeriod, "CARBONIFEROUS")
  expect_equal(byLatLong$maxAge, 313.5)
  bng <- LatLongToBNG(54.7733, -1.5762)
  expect_equal(Bedrock(easting = bng$easting, northing = bng$northing),
               byLatLong)
  centre <- GridRefToBNG("NZ 274 422")
  expect_equal(Bedrock(gridRef = "NZ 274 422"),
               Bedrock(easting = centre$easting, northing = centre$northing))
})

test_that("each map reports its own units", {
  # South Downs, Durham, Dogger Bank
  lat <- c(50.745, 54.7733, 55)
  long <- c(0.22, -1.5762, 2.5)
  bed <- Bedrock(lat, long)
  expect_equal(bed$rock[1], "CHALK")
  expect_match(bed$rock[2], "COAL")
  expect_true(is.na(bed$unit[3]))

  expect_equal(Superficial(gridRef = "NZ 274 422")$unit, "TILL")
  expect_true(all(is.na(Superficial(lat, long)$maxAge)))

  sea <- Seabed(lat, long)
  expect_true(is.na(sea$unit[2]))
  expect_equal(sea$folk[3], "S")
  expect_named(sea, c("easting", "northing", "lex", "unit", "rcs", "rock",
                      "maxPeriod", "minPeriod", "maxAge", "minAge", "folk"))
})

test_that("points off the map give NA", {
  res <- Bedrock(easting = c(NA, 427400, 427400, -2e6),
                 northing = c(542200, NA, 9e6, 542200))
  expect_equal(nrow(res), 4)
  expect_true(all(is.na(res$unit)))

  # Southernmost band of the index
  index <- .Layer("seabed")
  south <- (index$bandMin + 0.5) * .bandHeight * index$resolution
  expect_equal(nrow(Seabed(easting = 400000, northing = south)), 1)

  expect_equal(nrow(Bedrock(numeric(0), numeric(0))), 0)
})

test_that("swapped latitude and longitude prompt a warning", {
  expect_warning(res <- Bedrock(-1.5762, 54.7733), "right way round")
  expect_true(is.na(res$unit))
})

test_that("points must be given one way", {
  expect_error(Bedrock(), "Give `lat` and `long`")
  expect_error(Bedrock(54), "Give `lat` and `long`")
  expect_error(Bedrock(54, -1, easting = 1), "Give `lat` and `long`")
  expect_error(Bedrock(easting = 427400), "Give `lat` and `long`")
  expect_error(Bedrock(easting = "427400", northing = 542200), "numeric")
})

test_that("BGSSources() describes each map", {
  src <- BGSSources()
  expect_equal(src[["function"]], c("Bedrock", "Superficial", "Seabed"))
  expect_equal(src$scale, c(625000, 625000, 250000))
  expect_true(all(src$licence == "Open Government Licence v3.0"))
})
