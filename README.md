# BGS

What rock is under my feet? 'BGS' reports the bedrock, superficial deposits
and seabed sediments that the British Geological Survey maps at any point in
the UK or on its continental shelf. Give it latitude and longitude, British
National Grid coordinates, or an Ordnance Survey grid reference. The maps
are built in, so you need no GIS software and no internet connection.

## Install

```r
install.packages(
  "BGS",
  repos = c("https://geol-skills.github.io/packages", getOption("repos"))
)
```

In [webR](https://docs.r-wasm.org/webr/latest/):

```r
webr::install("BGS", repos = "https://geol-skills.github.io/packages")
```

## Use

```r
library("BGS")

# Durham Cathedral
Bedrock(54.7733, -1.5762)
Superficial(easting = 427400, northing = 542200)
Bedrock(gridRef = "NZ 274 422")

# Dogger Bank
Seabed(55, 2.5)
```

Each function returns a data frame with one row per point: the mapped unit,
its rock type, and the geological periods and ages it spans.

## Data

| Function | Map | Scale |
|---|---|---|
| `Bedrock()` | BGS Geology 625k (DiGMapGB-625), bedrock | 1:625 000 |
| `Superficial()` | BGS Geology 625k (DiGMapGB-625), superficial deposits | 1:625 000 |
| `Seabed()` | BGS Seabed Sediments 250k (DiGSBS250k), Folk classes | 1:250 000 |

Contains British Geological Survey materials © UKRI, released under the
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).
`BGSSources()` gives each map's version. `data-raw/` rebuilds the data from
the BGS downloads.

This package is not affiliated with or endorsed by the British Geological
Survey.
