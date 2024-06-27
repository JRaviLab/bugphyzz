## Script to generate sysdata.Rda

library(taxizedb)
library(bugphyzz)
library(purrr)
library(dplyr)
library(magrittr)
library(stringr)

getParentRank <- function(x) {
  ranks <- taxizedb::taxid2rank(x, db = 'ncbi', verbose = FALSE)
  lowestRanks <- c(
    'biotype', 'isolate', 'serogroup', 'serotype', 'strain', 'subspecies'
  )
  dplyr::case_when(
    ranks %in% lowestRanks ~ 'species',
    ranks == 'species' ~ 'genus',
    ranks == 'genus' ~ 'family',
    TRUE ~ NA
  )
}

taxRanks <- c(
  "superkingdom", "phylum", "class", "order", "family", "genus",
  "species", "strain"
)

phys <- physiologies()

ncbiIds <- phys |>
  map( ~ pull(.x, NCBI_ID)) |>
  flatten_chr() |>
  unique() |>
  str_squish() |>
  str_to_lower() |>
  {\(y) y[y != 'unknown']}() |>
  {\(y) y[!is.na(y)]}() |>
  sort(decreasing = TRUE)

tim <- system.time({
  taxonomies <- taxizedb::classification(ncbiIds, db = "ncbi")
  lglVct <- !map_lgl(taxonomies, ~ all(is.na(.x)))
  taxonomies <- taxonomies[lglVct]
  ncbiIds <- ncbiIds[lglVct]
})
print(tim)

## Check names and taxid match
all(names(taxonomies) == map_chr(taxonomies, ~ as.character(tail(.x$id, 1))))

parentsRanks <- getParentRank(ncbiIds)
lglVct <- !is.na(parentsRanks)
ncbiIds <- ncbiIds[lglVct]
parentsRanks <- parentsRanks[lglVct]
taxonomies <- taxonomies[lglVct]

parentIds <- map2(taxonomies, parentsRanks,  ~{
  parentRank <- .x |>
    filter(rank %in% taxRanks) |>
    pull(rank) |>
    {\(y) y[-length(y)]}() |> ## Need to remove the current rank
    tail(1)
  parentId <- .x |>
    filter(rank %in% taxRanks) |>
    pull(id) |>
    {\(y) y[-length(y)]}() |> ## Need to remove the current rank
    tail(1)
  names(parentId) <- parentRank
  ifelse(names(parentId) == .y, parentId, NA)
})

lglVct <- !is.na(parentIds)
ncbiIds <- ncbiIds[lglVct]
parentIds <- parentIds[lglVct]

ranksParents <- data.frame(
  NCBI_ID = ncbiIds,
  # Taxon_name = taxizedb::taxid2name(ncbiIds, db = 'ncbi'),
  Rank = taxizedb::taxid2rank(ncbiIds, db = 'ncbi'),
  Parent_NCBI_ID = unlist(parentIds),
  Parent_name = taxizedb::taxid2name(unlist(parentIds), db = 'ncbi'),
  Parent_rank = taxizedb::taxid2rank(unlist(parentIds), db = 'ncbi')
)
rownames(ranksParents) <- NULL

# BacDive -----------------------------------------------------------------
bacdive <- bugphyzz:::.getBacDive() |>
  bugphyzz:::.reshapeBacDive()
bacdivePhysNames <- names(bacdive)

## Save data -------------------------------------------------------------
usethis::use_data(
  ranksParents,
  bacdivePhysNames,
  overwrite = TRUE, internal = TRUE
)
