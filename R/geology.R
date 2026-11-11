#' Geology at a point
#'
#' `Bedrock()`, `Superficial()` and `Seabed()` report the bedrock, superficial
#' deposits and seabed sediments mapped by the British Geological Survey at
#' each point.
#'
#' Locate points by latitude and longitude, by British National Grid
#' coordinates, or by Ordnance Survey grid references.
#'
#' Bedrock and superficial deposits come from the 1:625 000 geological map of
#' the UK (onshore only); seabed sediments, from the 1:250 000 map of the UK
#' continental shelf. At these scales a boundary may be hundreds of metres from
#' its true position, so treat results near a boundary with caution.
#' `NA` means that nothing is mapped at that point: the sea, for `Bedrock()`;
#' ground without mapped superficial deposits, for `Superficial()`; land, or
#' sea beyond the map, for `Seabed()`.
#'
#' @inheritParams LatLongToBNG
#' @param gridRef Character vector of Ordnance Survey grid references, such as
#' `"NZ 274 422"`, each read as the centre of its grid square.
#'
#' @return A data frame with one row per point, and columns:
#' - `easting`, `northing`: location on the British National Grid (m);
#' - `lex`, `unit`: code and name of the mapped unit in the
#'   [BGS Lexicon](https://webapps.bgs.ac.uk/lexicon/);
#' - `rcs`, `rock`: code and description of its rock type;
#' - `maxPeriod`, `minPeriod`: oldest and youngest geological period it spans;
#' - `maxAge`, `minAge`: oldest and youngest age it spans, in millions of
#'   years (`NA` for superficial deposits);
#' - `folk` (`Seabed()` only): sediment class of Folk (1954), such as `"gmS"`
#'   for gravelly muddy sand.
#'
#' @references
#' Folk, R. L. (1954). The distinction between grain size and mineral
#' composition in sedimentary-rock nomenclature.
#' _Journal of Geology_ 62, 344--359. \doi{10.1086/626171}
#'
#' @examples
#' # Durham Cathedral
#' Bedrock(54.7733, -1.5762)
#' Superficial(easting = 427400, northing = 542200)
#'
#' # Dogger Bank
#' Seabed(55, 2.5)[, c("unit", "rock", "folk")]
#' @seealso [BGSSources()] for the maps' versions and licence.
#' @export
Bedrock <- function(lat, long, easting, northing, gridRef) {
  .Lookup("bedrock", .Locate(lat, long, easting, northing, gridRef))
}

#' @rdname Bedrock
#' @export
Superficial <- function(lat, long, easting, northing, gridRef) {
  .Lookup("superficial", .Locate(lat, long, easting, northing, gridRef))
}

#' @rdname Bedrock
#' @export
Seabed <- function(lat, long, easting, northing, gridRef) {
  .Lookup("seabed", .Locate(lat, long, easting, northing, gridRef))
}

#' Map sources
#'
#' Lists the British Geological Survey maps that this package contains, so
#' that you can cite the version you used.
#'
#' Contains British Geological Survey materials © UKRI, released under the
#' [Open Government Licence v3.0](
#' https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
#'
#' @return A data frame with one row per map, giving the function that reads
#' it, its title, version, release date, scale, licence and web page.
#' @examples
#' BGSSources()
#' @export
BGSSources <- function() {
  layers <- c(Bedrock = "bedrock", Superficial = "superficial",
              Seabed = "seabed")
  sources <- lapply(layers, function(name) {
    src <- .Layer(name)[["source"]]
    data.frame(
      title = src[["title"]],
      version = src[["version"]],
      released = src[["released"]],
      scale = src[["scale"]],
      licence = src[["licence"]],
      url = src[["url"]]
    )
  })
  cbind(`function` = names(layers), do.call(rbind, unname(sources)))
}

.Locate <- function(lat, long, easting, northing, gridRef) {
  given <- c(
    lat = !missing(lat), long = !missing(long),
    easting = !missing(easting), northing = !missing(northing),
    gridRef = !missing(gridRef)
  )
  form <- if (identical(unname(given), c(TRUE, TRUE, FALSE, FALSE, FALSE))) {
    "latLong"
  } else if (identical(unname(given), c(FALSE, FALSE, TRUE, TRUE, FALSE))) {
    "bng"
  } else if (identical(unname(given), c(FALSE, FALSE, FALSE, FALSE, TRUE))) {
    "gridRef"
  } else {
    stop("Give `lat` and `long`, or `easting` and `northing`, or `gridRef`.")
  }
  switch(form,
    latLong = {
      bng <- LatLongToBNG(lat, long)
      if (length(lat) && all(
        lat < 48 | lat > 64 | long < -14 | long > 5, na.rm = TRUE
      )) {
        warning(
          "No point lies near the UK: are `lat` and `long` the right way ",
          "round?"
        )
      }
      bng
    },
    bng = {
      .CheckPair(easting, northing, "easting", "northing")
      data.frame(easting = easting, northing = northing)
    },
    gridRef = GridRefToBNG(gridRef)
  )
}

.Lookup <- function(layer, bng) {
  index <- .Layer(layer)
  res <- index[["resolution"]]
  unit <- .Find(bng[["easting"]] / res, bng[["northing"]] / res, index)
  units <- index[["units"]][unit, , drop = FALSE]
  rownames(units) <- NULL
  cbind(bng, units)
}

.layers <- new.env(parent = emptyenv())

.Layer <- function(name) {
  if (is.null(.layers[[name]])) {
    file <- system.file("extdata", paste0(name, ".rds"), package = "BGS")
    .layers[[name]] <- .Index(readRDS(file))
  }
  .layers[[name]]
}

# Rows of the band index, in units of the layer's resolution
.bandHeight <- 400L

# Index each edge by every horizontal band its y-range touches, so that a ray
# cast east from a point need only test the edges in that point's band.
.Index <- function(layer) {
  x <- cumsum(layer[["dx"]])
  y <- cumsum(layer[["dy"]])
  ringEnd <- cumsum(layer[["ringLength"]])
  nxt <- seq_along(x) + 1L
  nxt[ringEnd] <- ringEnd - layer[["ringLength"]] + 1L
  lo <- pmin(y, y[nxt]) %/% .bandHeight
  hi <- pmax(y, y[nxt]) %/% .bandHeight
  bandMin <- min(lo)
  copies <- hi - lo + 1L
  band <- sequence(copies, from = lo - bandMin + 1L)
  edges <- rep.int(seq_along(x), copies)[order(band)]
  list(
    x = x, y = y, nxt = nxt,
    unit = rep.int(layer[["ringUnit"]], layer[["ringLength"]]),
    edges = edges,
    bandEnd = cumsum(tabulate(band)),
    bandMin = bandMin,
    units = layer[["units"]],
    resolution = layer[["resolution"]],
    source = layer[["source"]]
  )
}

# Index of the unit containing each point (px, py), by the even-odd rule
.Find <- function(px, py, index) {
  nUnit <- nrow(index[["units"]])
  bandEnd <- index[["bandEnd"]]
  x <- index[["x"]]
  y <- index[["y"]]
  band <- py %/% .bandHeight - index[["bandMin"]] + 1L
  vapply(seq_along(px), function(i) {
    b <- band[i]
    if (is.na(b) || is.na(px[i]) || b < 1L || b > length(bandEnd)) {
      return(NA_integer_)
    }
    from <- if (b == 1L) 1L else bandEnd[b - 1L] + 1L
    edge <- index[["edges"]][seq.int(from, length.out = bandEnd[b] - from + 1L)]
    to <- index[["nxt"]][edge]
    y1 <- y[edge]
    y2 <- y[to]
    spans <- (y1 > py[i]) != (y2 > py[i])
    edge <- edge[spans]
    to <- to[spans]
    y1 <- y1[spans]
    x1 <- x[edge]
    crossing <- x1 + (py[i] - y1) * (x[to] - x1) / (y2[spans] - y1)
    odd <- which(tabulate(index[["unit"]][edge[crossing > px[i]]], nUnit) %% 2L
                 == 1L)
    if (length(odd)) odd[1] else NA_integer_
  }, integer(1))
}
