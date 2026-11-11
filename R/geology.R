#' Geology at a point
#'
#' `Bedrock()` and `Superficial()` report the bedrock and the superficial
#' deposits mapped by the British Geological Survey at each point, on land or
#' under the sea.
#'
#' Locate points by latitude and longitude, by British National Grid
#' coordinates, or by Ordnance Survey grid references.
#'
#' On land, both functions read the 1:625 000 geological map of the UK; at
#' sea (anywhere off that map), the 1:250 000 maps of the UK continental shelf,
#' whose superficial deposits are seabed sediments.
#' At these scales a boundary may lie hundreds of metres from its true
#' position, so treat results near a boundary, or the coast, with caution.
#' `NA` means that nothing is mapped at that point: beyond the shelf maps, or,
#' for `Superficial()`, ground without mapped superficial deposits.
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
#'   years (`NA` for superficial deposits on land);
#' - `folk` (`Superficial()` only): sediment class of Folk (1954), such as
#'   `"gmS"` for gravelly muddy sand, for seabed sediments;
#' - `scale`: scale of the map that gave the result.
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
#' Bedrock(gridRef = "NZ 274 422")
#'
#' # Dogger Bank
#' Superficial(55, 2.5)[, c("unit", "rock", "folk")]
#' @seealso [BGSSources()] for the maps' versions and licence.
#' @export
Bedrock <- function(lat, long, easting, northing, gridRef) {
  .Geology(.Locate(lat, long, easting, northing, gridRef), "bedrock")
}

#' @rdname Bedrock
#' @export
Superficial <- function(lat, long, easting, northing, gridRef) {
  .Geology(.Locate(lat, long, easting, northing, gridRef), "superficial")
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
#' it, where it applies, and the map's title, version, release date, scale,
#' licence and web page.
#' @examples
#' BGSSources()
#' @export
BGSSources <- function() {
  maps <- expand.grid(
    covers = c("land", "sea"), layer = c("bedrock", "superficial"),
    stringsAsFactors = FALSE
  )
  sources <- lapply(paste0(maps$layer, "-", maps$covers), function(name) {
    as.data.frame(.Map(name)[["source"]])
  })
  cbind(
    `function` = c(bedrock = "Bedrock", superficial = "Superficial")[
      maps$layer
    ],
    covers = maps$covers,
    do.call(rbind, sources),
    row.names = NULL
  )
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
        lat < 48 | lat > 64 | long < -16 | long > 5, na.rm = TRUE
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

# Land is wherever the onshore bedrock map has a unit; elsewhere, use the
# offshore map.
.Geology <- function(bng, layer) {
  onshore <- .Find(bng, "bedrock-land")
  atSea <- is.na(onshore)
  landName <- paste0(layer, "-land")
  seaName <- paste0(layer, "-sea")
  land <- .Map(landName)
  sea <- .Map(seaName)
  landUnit <- if (layer == "bedrock") onshore else .Find(bng, landName)
  seaUnit <- .Find(bng[atSea, , drop = FALSE], seaName)
  row <- landUnit
  row[atSea] <- nrow(land[["units"]]) + seaUnit
  scale <- ifelse(atSea, sea[["source"]][["scale"]],
                  land[["source"]][["scale"]])
  scale[is.na(row)] <- NA
  units <- rbind(land[["units"]], sea[["units"]])[row, , drop = FALSE]
  rownames(units) <- NULL
  cbind(bng, units, scale = scale)
}

.maps <- new.env(parent = emptyenv())

.Map <- function(name) {
  if (is.null(.maps[[name]])) {
    file <- system.file("extdata", paste0(name, ".rds"), package = "BGS")
    .maps[[name]] <- .Index(readRDS(file))
  }
  .maps[[name]]
}

# Rows of the band index, in units of the map's resolution
.bandHeight <- 200L

# Index each edge by every horizontal band its y-range touches, so that a ray
# cast east from a point need only test the edges in that point's band.
.Index <- function(map) {
  x <- cumsum(map[["dx"]])
  y <- cumsum(map[["dy"]])
  ringEnd <- cumsum(map[["ringLength"]])
  nxt <- seq_along(x) + 1L
  nxt[ringEnd] <- ringEnd - map[["ringLength"]] + 1L
  lo <- pmin(y, y[nxt]) %/% .bandHeight
  hi <- pmax(y, y[nxt]) %/% .bandHeight
  bandMin <- min(lo)
  copies <- hi - lo + 1L
  band <- sequence(copies, from = lo - bandMin + 1L)
  unit <- rep.int(map[["ringUnit"]], map[["ringLength"]])
  map[c("dx", "dy", "ringLength", "ringUnit")] <- NULL
  c(map, list(
    x = x, y = y, nxt = nxt, unit = unit,
    edges = rep.int(seq_along(x), copies)[order(band)],
    bandEnd = cumsum(tabulate(band)),
    bandMin = bandMin
  ))
}

# Row of the map's `units` containing each point, by the even-odd rule
.Find <- function(bng, name) {
  map <- .Map(name)
  px <- bng[["easting"]] / map[["resolution"]]
  py <- bng[["northing"]] / map[["resolution"]]
  nUnit <- nrow(map[["units"]])
  bandEnd <- map[["bandEnd"]]
  x <- map[["x"]]
  y <- map[["y"]]
  band <- py %/% .bandHeight - map[["bandMin"]] + 1L
  vapply(seq_along(px), function(i) {
    b <- band[i]
    if (is.na(b) || is.na(px[i]) || b < 1L || b > length(bandEnd)) {
      return(NA_integer_)
    }
    from <- if (b == 1L) 1L else bandEnd[b - 1L] + 1L
    edge <- map[["edges"]][seq.int(from, length.out = bandEnd[b] - from + 1L)]
    to <- map[["nxt"]][edge]
    y1 <- y[edge]
    y2 <- y[to]
    spans <- (y1 > py[i]) != (y2 > py[i])
    edge <- edge[spans]
    to <- to[spans]
    y1 <- y1[spans]
    x1 <- x[edge]
    crossing <- x1 + (py[i] - y1) * (x[to] - x1) / (y2[spans] - y1)
    odd <- which(tabulate(map[["unit"]][edge[crossing > px[i]]], nUnit) %% 2L
                 == 1L)
    if (length(odd)) odd[1] else NA_integer_
  }, integer(1))
}
