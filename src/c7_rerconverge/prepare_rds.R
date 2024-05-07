# if (!require("devtools", quietly = TRUE)) install.packages("devtools")
# install_github("nclark-lab/RERconverge")

library(conflicted)
library(tidyverse)
library(RERconverge)


# Read trees -------------------------------------------------------------------
trees_nogap = RERconverge::readTrees(
  "out/genetree.pep.nogap.trees",
  # max.read = 10000,
  # minTreesAll = 19,
  # minSpecs = NULL
)
base::saveRDS(trees_nogap, file="rstat/rds/trees.nogap.rds")
ape::write.tree(trees_nogap$masterTree, file="out/RERconverge_master_tree.nogap.nwk")

trees_autom = RERconverge::readTrees(
  "out/genetree.pep.autom.trees",
  # max.read = 10000,
  # minTreesAll = 19,
  # minSpecs = NULL
)
base::saveRDS(trees_autom, file="rstat/rds/trees.autom.rds")
ape::write.tree(trees_autom$masterTree, file="out/RERconverge_master_tree.autom.nwk")


# Calculate relative evolutionary rates ----------------------------------------
lht = readr::read_tsv("out/species.tsv") |>
  dplyr::mutate(
    rownm = stringr::str_replace(`Organism Name`, " ", "_")
  ) |>  # space is not allowed in rownames?
  column_to_rownames(var = "rownm")
conflicts_prefer(stats::lowess)


## Mammals + Reptilia + Aves
lht_fore_aves = lht |>
  dplyr::filter(class  %in% c("Aves", "Reptilia", "Mammalia"))
RER_aves_nogap = RERconverge::getAllResiduals(
  trees_nogap,
  useSpecies = rownames(lht_fore_aves),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER_aves_nogap, file="rstat/rds/RER_fore_aves.nogap.rds")

RER_avse_autom = RERconverge::getAllResiduals(
  trees_autom,
  useSpecies = rownames(lht_fore_aves),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER, file="rstat/rds/RER_fore_aves.autom.rds")

## Mammals + Reptilia + Bats
lht_fore_bats = lht |>
  dplyr::filter(class  %in% c("Chiroptera", "Reptilia", "Mammalia"))
RER_bats_nogap = RERconverge::getAllResiduals(
  trees_nogap,
  useSpecies = rownames(lht_fore_bats),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER_bats_nogap, file="rstat/rds/RER_fore_bats.nogap.rds")

RER_bats_autom = RERconverge::getAllResiduals(
  trees_autom,
  useSpecies = rownames(lht_fore_bats),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER_bats_autom, file="rstat/rds/RER_fore_bats.autom.rds")


# ------------------------------------------------------------------------------

## All
RER_nogap = RERconverge::getAllResiduals(
  trees_nogap,
  useSpecies = rownames(lht),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER_nogap, file="rstat/rds/RER_all.nogap.rds")

RER_autom = RERconverge::getAllResiduals(
  trees_autom,
  useSpecies = rownames(lht),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER_autom, file="rstat/rds/RER_all.autom.rds")

## Mammals
lht_mammals = lht |>
  dplyr::filter(class %in% c("Mammalia", "Chiroptera"))
RER = RERconverge::getAllResiduals(
  trees,
  useSpecies = rownames(lht_mammals),
  transform = "sqrt",
  min.sp = 3,
  plot = FALSE
)
base::saveRDS(RER, file="rstat/rds/RER_mammals.rds")
