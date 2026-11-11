# Checks the package against terra (PROJ and GEOS). Run from the package root
# after build-layers.R, with the same download directory:
#   Rscript data-raw/validate.R [download-directory]

library("terra")
devtools::load_all(quiet = TRUE)

args <- commandArgs(TRUE)
dir <- if (length(args)) args[1] else file.path(tempdir(), "bgs-sources")
onshore <- file.path(dir, "625k_V5_Geology_UK_EPSG27700.gpkg")
offshore <- file.path(dir, "BGS_BedrockOffshore_250k_WGS84_v3.gpkg")
seabed <- file.path(dir, "BGS_250k_SeabedSediments_WGS84_v3_FOLK.gpkg")
set.seed(0)

# Coordinates: a grid over the UK, Ireland and the continental shelf
ll <- expand.grid(lat = seq(49, 61, by = 0.5), long = seq(-11, 3, by = 0.5))
ours <- LatLongToBNG(ll$lat, ll$long)
proj <- crds(project(vect(ll, geom = c("long", "lat"), crs = "EPSG:4326"),
                     "EPSG:27700"))
gap <- sqrt((ours$easting - proj[, 1]) ^ 2 + (ours$northing - proj[, 2]) ^ 2)
# Beyond its OSTN15 grid, PROJ falls back to a "ballpark" transformation with
# no datum shift, so compare only on mapped British land
land <- !is.na(.Find(ours, "bedrock-land"))
cat("LatLongToBNG vs PROJ on land: median", signif(median(gap[land]), 3),
    "m; max", signif(max(gap[land]), 3), "m\n")
back <- BNGToLatLong(ours$easting, ours$northing)
cat("Round trip: max error", signif(max(abs(c(
  back$lat - ll$lat, back$long - ll$long
))) * 111e3, 3), "m\n")

# Lookups: random points over each map, compared with terra::extract on the
# map as published (unrounded; the 250k maps in WGS 84, so that this also
# tests LatLongToBNG() across the shelf)
Compare <- function(v, name, fields, n = 1e5) {
  e <- ext(v)
  pts <- data.frame(x = runif(n, e$xmin, e$xmax), y = runif(n, e$ymin, e$ymax))
  bng <- if (is.lonlat(v)) {
    LatLongToBNG(pts$y, pts$x)
  } else {
    data.frame(easting = pts$x, northing = pts$y)
  }
  found <- .Map(name)$units[.Find(bng, name), ]
  mine <- do.call(paste, found[c("lex", "rcs")])
  mine[is.na(found$lex)] <- NA
  where <- vect(pts, geom = c("x", "y"), crs = crs(v))
  hit <- extract(v, where)
  hit <- hit[!duplicated(hit$id.y), ]
  theirs <- do.call(paste, lapply(fields, function(f) trimws(hit[[f]])))
  theirs[is.na(hit[[fields[1]]])] <- NA
  differ <- which(xor(is.na(mine), is.na(theirs)) |
                    (!is.na(mine) & !is.na(theirs) & mine != theirs))
  cat(name, ": ", length(differ), " of ", n, " points differ", sep = "")
  if (length(differ)) {
    near <- distance(where[differ], as.lines(v))
    cat("; furthest from a boundary:", signif(max(apply(near, 1, min)), 3),
        "m")
  }
  cat("\n")
}

Compare(vect(onshore, layer = "625k_V5_BEDROCK_Geology"), "bedrock-land",
        c("LEX", "RCS"))
Compare(vect(offshore, layer = vector_layers(offshore)[1]), "bedrock-sea",
        c("LEX", "RCS"))
Compare(vect(onshore, layer = "625k_V5_SUPERFICIAL_Geology"),
        "superficial-land", c("LEX", "ROCK"))
Compare(vect(seabed, layer = vector_layers(seabed)[1]), "superficial-sea",
        c("LEX", "RCS"))
