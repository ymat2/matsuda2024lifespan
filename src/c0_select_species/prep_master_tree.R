library(conflicted)
library(tidyverse)
library(phytools)


# asMasterTree -----------------------------------------------------------------
rename_tip = function(phy, old_name, new_name) {
  mpos = base::match(old_name, phy$tip.label)
  phy$tip.label[mpos] = new_name
  return(phy)
}

monoBranch = function(phy) {
  len = length(phy$edge.length)
  x = rep(c(1), length = len)
  phy$edge.length = x
  return(phy)
}

tre = ape::read.tree("data/species.nwk")
df_ncbi = readr::read_tsv("out/species.tsv")

tips = base::intersect(
  tre$tip.label,
  stringr::str_replace_all(df_ncbi$`Organism Name`, " ", "_")
)
tre = ape::keep.tip(tre, tips)

master_tree = monoBranch(tre)
ape::write.tree(master_tree, file = "out/mastertree.nwk")


# filter table and species -----------------------------------------------------
df_ncbi_on_tree = df_ncbi |>
  dplyr::filter(`Organism Name` %in% stringr::str_replace_all(tips, "_", " "))

readr::write_tsv(df_ncbi_on_tree, file = "out/species.tsv")
readr::write_lines(df_ncbi_on_tree$accession, "out/accession.list")
