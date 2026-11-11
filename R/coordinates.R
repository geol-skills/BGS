#' Convert between latitude and longitude, grid coordinates and grid references
#'
#' `LatLongToBNG()` converts latitude and longitude to British National Grid
#' eastings and northings; `BNGToLatLong()` converts back.
#' `GridRefToBNG()` reads Ordnance Survey grid references.
#'
#' Latitude and longitude are taken to be WGS 84, as used by GPS and web maps.
#' Conversion to the OSGB36 datum of the National Grid uses the Ordnance
#' Survey's seven-parameter Helmert transformation, accurate to about 5 m:
#' far finer than the maps in this package.
#'
#' @param lat,long Numeric vectors giving latitude and longitude in decimal
#' degrees, north and east positive.
#' @param easting,northing Numeric vectors giving British National Grid
#' coordinates in metres.
#'
#' @return `LatLongToBNG()` and `GridRefToBNG()` return a data frame with
#' columns `easting` and `northing`; `BNGToLatLong()`, a data frame with
#' columns `lat` and `long`.
#'
#' @references
#' Ordnance Survey (2020). _A guide to coordinate systems in Great Britain_.
#' Southampton: Ordnance Survey.
#'
#' @examples
#' # Durham Cathedral
#' LatLongToBNG(54.7733, -1.5762)
#' BNGToLatLong(427400, 542200)
#' GridRefToBNG(c("NZ 274 422", "NZ2742"))
#' @family coordinate conversions
#' @export
LatLongToBNG <- function(lat, long) {
  .CheckPair(lat, long, "lat", "long")
  if (any(abs(lat) > 90, na.rm = TRUE)) {
    stop("`lat` must lie between -90 and 90 degrees.")
  }
  if (any(abs(long) > 180, na.rm = TRUE)) {
    stop("`long` must lie between -180 and 180 degrees.")
  }
  xyz <- .ToCartesian(lat, long, .grs80)
  osgb <- .FromCartesian(.Helmert(xyz, .wgs84ToOsgb), .airy)
  .TransverseMercator(osgb[["lat"]], osgb[["long"]])
}

#' @rdname LatLongToBNG
#' @export
BNGToLatLong <- function(easting, northing) {
  .CheckPair(easting, northing, "easting", "northing")
  osgb <- .InverseTransverseMercator(easting, northing)
  xyz <- .ToCartesian(osgb[["lat"]], osgb[["long"]], .airy)
  .FromCartesian(.Helmert(xyz, -.wgs84ToOsgb), .grs80)
}

#' @rdname LatLongToBNG
#' @param gridRef Character vector of Ordnance Survey grid references: two
#' letters then an even number of digits, such as `"NZ 274 422"`; spaces and
#' case are ignored.
#' @param centre Logical; if `TRUE`, return the centre of each grid square,
#' rather than its south-west corner.
#' @export
GridRefToBNG <- function(gridRef, centre = TRUE) {
  ref <- toupper(gsub("[[:space:]]", "", as.character(gridRef)))
  bad <- !is.na(ref) & !grepl("^[HJNOST][A-HJ-Z]([0-9][0-9]){0,5}$", ref)
  if (any(bad)) {
    stop(
      "Can't read grid reference ",
      paste0("\"", gridRef[bad], "\"", collapse = ", "),
      ": expected two letters and an even number of digits, e.g. \"NZ 27 42\"."
    )
  }
  letter1 <- match(substr(ref, 1, 1), .gridLetters) - 1
  letter2 <- match(substr(ref, 2, 2), .gridLetters) - 1
  digits <- substring(ref, 3)
  half <- nchar(digits) / 2
  size <- 10 ^ (5 - half)
  eDigits <- as.numeric(substr(digits, 1, half))
  nDigits <- as.numeric(substring(digits, half + 1))
  eDigits[!is.na(ref) & half == 0] <- 0
  nDigits[!is.na(ref) & half == 0] <- 0
  offset <- if (isTRUE(centre)) size / 2 else 0
  data.frame(
    easting = ((letter1 - 2) %% 5 * 5 + letter2 %% 5) * 1e5 +
      eDigits * size + offset,
    northing = (19 - letter1 %/% 5 * 5 - letter2 %/% 5) * 1e5 +
      nDigits * size + offset
  )
}

.gridLetters <- LETTERS[LETTERS != "I"]

.CheckPair <- function(a, b, aName, bName) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("`", aName, "` and `", bName, "` must be numeric.")
  }
  if (length(a) != length(b)) {
    stop("`", aName, "` and `", bName, "` must have the same length.")
  }
}

.grs80 <- c(a = 6378137, b = 6356752.314140)
.airy <- c(a = 6377563.396, b = 6356256.909)

# tx, ty, tz (m); s (ppm); rx, ry, rz (arcseconds).  Negated, these give the
# reverse transformation to the same accuracy.
.wgs84ToOsgb <- c(
  tx = -446.448, ty = 125.157, tz = -542.060,
  s = 20.4894,
  rx = -0.1502, ry = -0.2470, rz = -0.8421
)

.Radians <- function(degrees) degrees * pi / 180

.ToCartesian <- function(lat, long, ellipsoid) {
  a <- ellipsoid[["a"]]
  e2 <- 1 - (ellipsoid[["b"]] / a) ^ 2
  phi <- .Radians(lat)
  lambda <- .Radians(long)
  nu <- a / sqrt(1 - e2 * sin(phi) ^ 2)
  list(
    x = nu * cos(phi) * cos(lambda),
    y = nu * cos(phi) * sin(lambda),
    z = (1 - e2) * nu * sin(phi)
  )
}

.Helmert <- function(xyz, p) {
  s1 <- 1 + p[["s"]] * 1e-6
  rx <- .Radians(p[["rx"]] / 3600)
  ry <- .Radians(p[["ry"]] / 3600)
  rz <- .Radians(p[["rz"]] / 3600)
  list(
    x = p[["tx"]] + s1 * xyz$x - rz * xyz$y + ry * xyz$z,
    y = p[["ty"]] + rz * xyz$x + s1 * xyz$y - rx * xyz$z,
    z = p[["tz"]] - ry * xyz$x + rx * xyz$y + s1 * xyz$z
  )
}

.FromCartesian <- function(xyz, ellipsoid) {
  a <- ellipsoid[["a"]]
  e2 <- 1 - (ellipsoid[["b"]] / a) ^ 2
  p <- sqrt(xyz$x ^ 2 + xyz$y ^ 2)
  phi <- atan2(xyz$z, p * (1 - e2))
  repeat {
    nu <- a / sqrt(1 - e2 * sin(phi) ^ 2)
    nextPhi <- atan2(xyz$z + e2 * nu * sin(phi), p)
    converged <- all(abs(nextPhi - phi) < 1e-12, na.rm = TRUE)
    phi <- nextPhi
    if (converged) break
  }
  data.frame(lat = phi * 180 / pi, long = atan2(xyz$y, xyz$x) * 180 / pi)
}

# National Grid projection constants
.f0 <- 0.9996012717
.phi0 <- .Radians(49)
.lambda0 <- .Radians(-2)
.e0 <- 400000
.n0 <- -100000

.Meridional <- function(phi) {
  a <- .airy[["a"]]
  b <- .airy[["b"]]
  n <- (a - b) / (a + b)
  dPhi <- phi - .phi0
  sPhi <- phi + .phi0
  b * .f0 * (
    (1 + n + 5 / 4 * n ^ 2 + 5 / 4 * n ^ 3) * dPhi -
      (3 * n + 3 * n ^ 2 + 21 / 8 * n ^ 3) * sin(dPhi) * cos(sPhi) +
      (15 / 8 * n ^ 2 + 15 / 8 * n ^ 3) * sin(2 * dPhi) * cos(2 * sPhi) -
      35 / 24 * n ^ 3 * sin(3 * dPhi) * cos(3 * sPhi)
  )
}

.Radii <- function(phi) {
  a <- .airy[["a"]]
  e2 <- 1 - (.airy[["b"]] / a) ^ 2
  w <- 1 - e2 * sin(phi) ^ 2
  nu <- a * .f0 / sqrt(w)
  rho <- a * .f0 * (1 - e2) / w ^ 1.5
  list(nu = nu, rho = rho, eta2 = nu / rho - 1)
}

# Latitude and longitude on the OSGB36 datum to National Grid
.TransverseMercator <- function(lat, long) {
  phi <- .Radians(lat)
  r <- .Radii(phi)
  nu <- r$nu
  s <- sin(phi)
  c <- cos(phi)
  t2 <- tan(phi) ^ 2
  dl <- .Radians(long) - .lambda0
  northing <- .Meridional(phi) + .n0 +
    nu / 2 * s * c * dl ^ 2 +
    nu / 24 * s * c ^ 3 * (5 - t2 + 9 * r$eta2) * dl ^ 4 +
    nu / 720 * s * c ^ 5 * (61 - 58 * t2 + t2 ^ 2) * dl ^ 6
  easting <- .e0 + nu * c * dl +
    nu / 6 * c ^ 3 * (nu / r$rho - t2) * dl ^ 3 +
    nu / 120 * c ^ 5 * (5 - 18 * t2 + t2 ^ 2 + 14 * r$eta2 -
                          58 * t2 * r$eta2) * dl ^ 5
  data.frame(easting = easting, northing = northing)
}

.InverseTransverseMercator <- function(easting, northing) {
  aF0 <- .airy[["a"]] * .f0
  phi <- (northing - .n0) / aF0 + .phi0
  repeat {
    residual <- northing - .n0 - .Meridional(phi)
    phi <- phi + residual / aF0
    if (all(abs(residual) < 1e-5, na.rm = TRUE)) break
  }
  r <- .Radii(phi)
  nu <- r$nu
  rho <- r$rho
  t <- tan(phi)
  t2 <- t ^ 2
  secPhi <- 1 / cos(phi)
  de <- easting - .e0
  xiiA <- secPhi / (5040 * nu ^ 7) *
    (61 + 662 * t2 + 1320 * t2 ^ 2 + 720 * t2 ^ 3)
  lat <- phi -
    t / (2 * rho * nu) * de ^ 2 +
    t / (24 * rho * nu ^ 3) * (5 + 3 * t2 + r$eta2 - 9 * t2 * r$eta2) * de ^ 4 -
    t / (720 * rho * nu ^ 5) * (61 + 90 * t2 + 45 * t2 ^ 2) * de ^ 6
  long <- .lambda0 +
    secPhi / nu * de -
    secPhi / (6 * nu ^ 3) * (nu / rho + 2 * t2) * de ^ 3 +
    secPhi / (120 * nu ^ 5) * (5 + 28 * t2 + 24 * t2 ^ 2) * de ^ 5 -
    xiiA * de ^ 7
  data.frame(lat = lat * 180 / pi, long = long * 180 / pi)
}
