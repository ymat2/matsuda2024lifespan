# devtools::install_github("nclark-lab/RERconverge")

library(conflicted)
library(tidyverse)
library(RERconverge)


# Reading RDS ------------------------------------------------------------------

trees_nogap = readr::read_rds("data/trees.nogap.rds")
trees_autom= readr::read_rds("data/trees.autom.rds")
RER_nogap = readr::read_rds("data/RER_all.nogap.rds")
RER_autom = readr::read_rds("data/RER_all.autom.rds")


# Binary trait analysis --------------------------------------------------------

lht = readr::read_tsv("out/species.tsv") |>
  dplyr::mutate(rownm = stringr::str_replace(`Organism Name`, " ", "_")) |>
  tibble::column_to_rownames(var = "rownm")


### nogap ----------------------------------------------------------------------

foreground = rownames(dplyr::filter(lht, class %in% c("Aves", "Chiroptera")))
volant = foreground2Tree(
  foreground,
  trees_nogap,
  clade="all",  # alternatively [all, terminal, ancestral]
  useSpecies=rownames(lht)
)
phenv=tree2Paths(volant, trees_nogap)

res_nogap=RERconverge::correlateWithBinaryPhenotype(
  RER_nogap,
  phenv,
  min.sp=10,
  min.pos=10,
  weighted="auto"
) |>
  tibble::rownames_to_column("symbol") |>
  dplyr::mutate(symbol = str_replace(symbol, ".pep.aln.*", "")) |>
  dplyr::arrange(P) |>
  tidyr::drop_na()

readr::write_tsv(res_nogap, "out/binary.nogap.tsv")


### automated1 -----------------------------------------------------------------

foreground = rownames(dplyr::filter(lht, class %in% c("Aves", "Chiroptera")))
volant = foreground2Tree(
  foreground,
  trees_autom,
  clade="all",  # alternatively [all, terminal, ancestral]
  useSpecies=rownames(lht)
)
phenv=tree2Paths(volant, trees_autom)

res_autom=RERconverge::correlateWithBinaryPhenotype(
  RER_autom,
  phenv,
  min.sp=10,
  min.pos=10,
  weighted="auto"
) |>
  tibble::rownames_to_column("symbol") |>
  dplyr::mutate(symbol = str_replace(symbol, ".pep.aln.*", "")) |>
  dplyr::arrange(P) |>
  tidyr::drop_na()

readr::write_tsv(res_autom, "out/binary.autom.tsv")
