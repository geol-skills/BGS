test_that("lookups agree whichever way a point is given", {
  # Durham Cathedral
  byLatLong <- Bedrock(54.7733, -1.5762)
  expect_equal(byLatLong$lex, "PSMCM")
  expect_equal(byLatLong$maxPeriod, "CARBONIFEROUS")
  expect_equal(byLatLong$maxAge, 313.5)
  expect_equal(byLatLong$scale, 625000)
  bng <- LatLongToBNG(54.7733, -1.5762)
  expect_equal(Bedrock(easting = bng$easting, northing = bng$northing),
               byLatLong)
  centre <- GridRefToBNG("NZ 274 422")
  expect_equal(Bedrock(gridRef = "NZ 274 422"),
               Bedrock(easting = centre$easting, northing = centre$northing))
})

test_that("land and sea read different maps", {
  # South Downs, Durham, Dogger Bank, mid-Atlantic
  lat <- c(50.745, 54.7733, 55, 50)
  long <- c(0.22, -1.5762, 2.5, -25)
  bed <- Bedrock(lat, long)
  expect_equal(bed$rock[1], "CHALK")
  expect_match(bed$rock[2], "COAL")
  expect_false(is.na(bed$unit[3]))
  expect_equal(bed$scale, c(625000, 625000, 250000, NA))
  expect_true(is.na(bed$unit[4]))
  expect_false("folk" %in% names(bed))

  sup <- Superficial(lat, long)
  expect_true(is.na(sup$unit[1]))
  expect_equal(sup$scale, c(NA, 625000, 250000, NA))
  expect_true(is.na(sup$folk[2]))
  expect_equal(sup$folk[3], "S")
  expect_true(is.na(sup$maxAge[2]))
  expect_equal(sup$maxAge[3], 0.01)
  expect_named(sup, c("easting", "northing", "lex", "unit", "rcs", "rock",
                      "maxPeriod", "minPeriod", "maxAge", "minAge", "folk",
                      "scale"))

  expect_equal(Superficial(gridRef = "NZ 274 422")$unit, "TILL")
})

test_that("points off the map give NA", {
  res <- Bedrock(easting = c(NA, 427400, 427400, -2e6),
                 northing = c(542200, NA, 9e6, 542200))
  expect_equal(nrow(res), 4)
  expect_true(all(is.na(res$unit)))

  # Southernmost band of the index
  map <- .Map("superficial-sea")
  south <- (map$bandMin + 0.5) * .bandHeight * map$resolution
  expect_equal(length(.Find(data.frame(easting = 4e5, northing = south),
                            "superficial-sea")), 1)

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
  expect_equal(src[["function"]],
               c("Bedrock", "Bedrock", "Superficial", "Superficial"))
  expect_equal(src$covers, c("land", "sea", "land", "sea"))
  expect_equal(src$scale, c(625000, 250000, 625000, 250000))
  expect_true(all(src$licence == "Open Government Licence v3.0"))
})
