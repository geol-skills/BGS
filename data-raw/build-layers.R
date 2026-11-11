# Builds inst/extdata/{bedrock,superficial}-{land,sea}.rds from the BGS
# GeoPackages. Run from the package root:
#   Rscript data-raw/build-layers.R [download-directory]
# Needs 'terra' (not a package dependency).
#
# Each map is a list:
#   units      data frame, one row per mapped unit, in the columns returned to
#              users (see R/geology.R); the land and sea maps of each function
#              share columns;
#   ringUnit   integer: the unit (row of `units`) each ring belongs to;
#   ringLength integer: vertices per ring, closing vertex dropped;
#   dx, dy     integer: vertex coordinates on the British National Grid, in
#              units of `resolution` metres, delta-encoded across all rings;
#   resolution metres per coordinate unit;
#   source     provenance, reported by BGSSources().
#
# Polygons of one unit never overlap, so a point lies in that unit exactly when
# a ray from it crosses the unit's rings an odd number of times: rings need not
# be grouped into polygons.

library("terra")
# The package's own conversion, so that stored and queried coordinates match:
# beyond the OSTN15 grid, PROJ's WGS 84 to OSGB36 omits the datum shift
source(file.path("R", "coordinates.R"))

resolution <- 10 # metres; 0.04 mm at 1:250 000

sources <- list(
  onshore = list(
    page = "https://www.bgs.ac.uk/download/bgs-geology-625k-bedrock/",
    zip = "BGS_Geology_625k_bedrock_gpkg.zip",
    gpkg = "625k_V5_Geology_UK_EPSG27700.gpkg",
    md5 = "559681730ddf9355884b39416210ee19"
  ),
  offshoreBedrock = list(
    page = paste0(
      "https://www.bgs.ac.uk/download/",
      "offshore-bedrock-250k-gis-geopackage/"
    ),
    zip = "offshore-bedrock-250k-geopackage.zip",
    gpkg = "BGS_BedrockOffshore_250k_WGS84_v3.gpkg",
    md5 = "e57994eb8b46129f01181fb93660943b"
  ),
  seabed = list(
    page = paste0(
      "https://www.bgs.ac.uk/download/",
      "offshore-seabed-sediments-250k-gis-geopackage/"
    ),
    zip = "bgs-offshore-sbs-250k-geopackage.zip",
    gpkg = "BGS_250k_SeabedSediments_WGS84_v3_FOLK.gpkg",
    md5 = "d058316b53dccb4671e2a449a3ff86ee"
  )
)

# BGS serves each file through a download page whose link carries a token
Fetch <- function(src, dir) {
  gpkg <- file.path(dir, src[["gpkg"]])
  if (!file.exists(gpkg)) {
    zip <- file.path(dir, src[["zip"]])
    if (!file.exists(zip)) {
      page <- readLines(src[["page"]], warn = FALSE)
      pattern <- paste0(
        "https://www\\.bgs\\.ac\\.uk/download/[^\"' ]*",
        "\\?wpdmdl=[0-9]+&(amp;)?refresh=[0-9a-f]+"
      )
      link <- regmatches(page, regexpr(pattern, page))[1]
      if (is.na(link)) stop("No download link on ", src[["page"]])
      download.file(gsub("&amp;", "&", link, fixed = TRUE), zip, mode = "wb")
    }
    unzip(zip, exdir = dir)
  }
  md5 <- unname(tools::md5sum(gpkg))
  if (md5 != src[["md5"]]) {
    warning(
      basename(gpkg), " has changed (MD5 ", md5, "): check the attributes ",
      "still match, then update `md5` in data-raw/build-layers.R"
    )
  }
  gpkg
}

Blank <- function(x) {
  x <- trimws(x)
  x[x %in% c("", "-", "UNKNOWN")] <- NA_character_
  x
}

Ma <- function(years) as.numeric(years) / 1e6

# Latest release date as YYYY-MM-DD; the 625k map writes "11- 2-2008"
Released <- function(x) {
  if (!inherits(x, c("Date", "POSIXt"))) {
    x <- as.Date(gsub(" ", "", x, fixed = TRUE), "%d-%m-%Y")
  }
  as.character(max(as.Date(x)))
}

# The 250k maps are in WGS 84, with edges straight in latitude and longitude
# (as GIS draws them) but some tens of kilometres long: add vertices so that
# they stay within a metre of that line once projected
ReadLatLong <- function(file) {
  densify(vect(file, layer = vector_layers(file)[1]), interval = 0.005,
          flat = TRUE)
}

Encode <- function(v, units, source) {
  key <- do.call(paste, units)
  unitTable <- units[!duplicated(key), , drop = FALSE]
  rownames(unitTable) <- NULL
  polyUnit <- match(key, key[!duplicated(key)])

  g <- geom(v)
  if (is.lonlat(v)) {
    bng <- LatLongToBNG(g[, "y"], g[, "x"])
    g[, "x"] <- bng[["easting"]]
    g[, "y"] <- bng[["northing"]]
  }
  x <- as.integer(round(g[, "x"] / resolution))
  y <- as.integer(round(g[, "y"] / resolution))
  newRing <- c(TRUE, diff(g[, "geom"]) != 0 | diff(g[, "part"]) != 0 |
                 diff(g[, "hole"]) != 0)
  ring <- cumsum(newRing)
  ringFirst <- which(newRing)[ring]

  # Drop each ring's closing vertex, and vertices that rounding has merged
  closing <- c(newRing[-1], TRUE) & x == x[ringFirst] & y == y[ringFirst]
  repeated <- c(FALSE, !newRing[-1] & diff(x) == 0 & diff(y) == 0)
  keep <- !closing & !repeated
  ringLength <- tabulate(ring[keep], nbins = max(ring))
  ringGeom <- g[newRing, "geom"]
  usable <- ringLength >= 3
  keep <- keep & usable[ring]

  x <- x[keep]
  y <- y[keep]
  list(
    units = unitTable,
    ringUnit = polyUnit[ringGeom[usable]],
    ringLength = ringLength[usable],
    dx = c(x[1], diff(x)),
    dy = c(y[1], diff(y)),
    resolution = resolution,
    source = c(source, licence = "Open Government Licence v3.0")
  )
}

Save <- function(map, name) {
  saveRDS(map, file.path(out, paste0(name, ".rds")), compress = "xz")
}

args <- commandArgs(TRUE)
dir <- if (length(args)) args[1] else file.path(tempdir(), "bgs-sources")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
files <- lapply(sources, Fetch, dir = dir)
out <- file.path("inst", "extdata")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
unlink(list.files(out, "\\.rds$", full.names = TRUE))
url625 <- "https://www.bgs.ac.uk/datasets/bgs-geology-625k-digmapgb/"

bed <- vect(files[["onshore"]], layer = "625k_V5_BEDROCK_Geology")
Save(Encode(bed, data.frame(
  lex = Blank(bed$LEX),
  unit = Blank(bed$LEX_D),
  rcs = Blank(bed$RCS),
  rock = Blank(bed$RCS_D),
  maxPeriod = Blank(bed$MAX_PERIOD),
  minPeriod = Blank(bed$MIN_PERIOD),
  maxAge = Ma(bed$MAX_TIME_Y),
  minAge = Ma(bed$MIN_TIME_Y)
), list(
  title = "BGS Geology 625k (DiGMapGB-625), bedrock",
  version = unique(bed$VERSION),
  released = Released(bed$RELEASED),
  scale = 625000,
  url = url625
)), "bedrock-land")

off <- ReadLatLong(files[["offshoreBedrock"]])
Save(Encode(off, data.frame(
  lex = Blank(off$LEX),
  unit = Blank(off$LEX_D),
  rcs = Blank(off$RCS),
  rock = Blank(off$RCS_D),
  maxPeriod = Blank(off$MAX_PERIOD),
  minPeriod = Blank(off$MIN_PERIOD),
  maxAge = Ma(off$MAX_TIME_Y),
  minAge = Ma(off$MIN_TIME_Y)
), list(
  title = "BGS Offshore Bedrock 250k (DiGRock250k)",
  version = unique(off$VERSION),
  released = Released(off$RELEASED),
  scale = 250000,
  url = "https://www.bgs.ac.uk/datasets/offshore-bedrock-250k/"
)), "bedrock-sea")

sup <- vect(files[["onshore"]], layer = "625k_V5_SUPERFICIAL_Geology")
Save(Encode(sup, data.frame(
  lex = Blank(sup$LEX),
  unit = Blank(sup$LEX_D),
  rcs = Blank(sup$ROCK),
  rock = Blank(sup$ROCK_D),
  maxPeriod = Blank(sup$MAX_SYSTEM),
  minPeriod = Blank(sup$MIN_SYSTEM),
  maxAge = NA_real_,
  minAge = NA_real_,
  folk = NA_character_
), list(
  title = "BGS Geology 625k (DiGMapGB-625), superficial deposits",
  version = unique(sup$VERSION),
  released = Released(sup$RELEASED),
  scale = 625000,
  url = url625
)), "superficial-land")

sea <- ReadLatLong(files[["seabed"]])
Save(Encode(sea, data.frame(
  lex = Blank(sea$LEX),
  unit = Blank(sea$LEX_D),
  rcs = Blank(sea$RCS),
  rock = Blank(sea$FOLK_D50),
  maxPeriod = Blank(sea$MAX_PERIOD),
  minPeriod = Blank(sea$MIN_PERIOD),
  maxAge = Ma(sea$MAX_TIME_Y),
  minAge = Ma(sea$MIN_TIME_Y),
  folk = Blank(sea$FOLK_S)
), list(
  title = "BGS Seabed Sediments 250k (DiGSBS250k), Folk classification",
  version = as.character(unique(sea$VERSION)),
  released = Released(sea$RELEASED),
  scale = 250000,
  url = "https://www.bgs.ac.uk/datasets/seabed-sediments-250k/"
)), "superficial-sea")

sizes <- file.size(list.files(out, full.names = TRUE)) / 1e6
print(setNames(round(sizes, 2), list.files(out)))
cat("Total:", round(sum(sizes), 2), "MB\n")
