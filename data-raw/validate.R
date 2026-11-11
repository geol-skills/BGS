# Checks the package against terra (PROJ and GEOS). Run from the package root
# after build-layers.R, with the same download directory:
#   Rscript data-raw/validate.R [download-directory]

library("terra")
devtools::load_all(quiet = TRUE)

args <- commandArgs(TRUE)
dir <- if (length(args)) args[1] else file.path(tempdir(), "bgs-sources")
onshore <- file.path(dir, "625k_V5_Geology_UK_EPSG27700.gpkg")
seabedFile <- file.path(dir, "BGS_250k_SeabedSediments_WGS84_v3_FOLK.gpkg")
set.seed(0)

# Coordinates: a grid over the UK, Ireland and the continental shelf
ll <- expand.grid(lat = seq(49, 61, by = 0.5), long = seq(-11, 3, by = 0.5))
ours <- LatLongToBNG(ll$lat, ll$long)
proj <- crds(project(vect(ll, geom = c("long", "lat"), crs = "EPSG:4326"),
                     "EPSG:27700"))
gap <- sqrt((ours$easting - proj[, 1]) ^ 2 + (ours$northing - proj[, 2]) ^ 2)
# Beyond its OSTN15 grid, PROJ falls back to a "ballpark" transformation with
# no datum shift, so compare only on mapped British land
land <- !is.na(Bedrock(easting = ours$easting, northing = ours$northing)$lex)
cat("LatLongToBNG vs PROJ on land: median", signif(median(gap[land]), 3),
    "m; max", signif(max(gap[land]), 3), "m\n")
back <- BNGToLatLong(ours$easting, ours$northing)
cat("Round trip: max error", signif(max(abs(c(
  back$lat - ll$lat, back$long - ll$long
))) * 111e3, 3), "m\n")

# Lookups: random points over each map, against terra::extract on the
# unrounded polygons
Compare <- function(v, Fun, fields, n = 1e5) {
  e <- ext(v)
  pts <- data.frame(
    easting = runif(n, e$xmin, e$xmax),
    northing = runif(n, e$ymin, e$ymax)
  )
  ours <- Fun(easting = pts$easting, northing = pts$northing)
  hit <- extract(v, vect(pts, geom = c("easting", "northing"),
                         crs = crs(v)))
  hit <- hit[!duplicated(hit$id.y), ]
  theirs <- do.call(paste, lapply(fields, function(f) trimws(hit[[f]])))
  theirs[is.na(hit[[fields[1]]])] <- NA
  mine <- do.call(paste, ours[c("lex", "rcs")])
  mine[is.na(ours$lex)] <- NA
  differ <- which(xor(is.na(mine), is.na(theirs)) |
                    (!is.na(mine) & !is.na(theirs) & mine != theirs))
  cat(deparse(substitute(Fun)), ": ", length(differ), " of ", n,
      " points differ", sep = "")
  if (length(differ)) {
    edges <- as.lines(v)
    near <- distance(vect(pts[differ, ], geom = c("easting", "northing"),
                          crs = crs(v)), edges)
    cat("; furthest from a boundary:", signif(max(apply(near, 1, min)), 3),
        "m")
  }
  cat("\n")
}

Compare(vect(onshore, layer = "625k_V5_BEDROCK_Geology"), Bedrock,
        c("LEX", "RCS"))
Compare(vect(onshore, layer = "625k_V5_SUPERFICIAL_Geology"), Superficial,
        c("LEX", "ROCK"))
Compare(project(vect(seabedFile), "EPSG:27700"), Seabed, c("LEX", "RCS"))
