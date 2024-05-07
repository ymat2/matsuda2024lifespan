library(conflicted)
library(tidyverse)


# Linear regression of Life history traits -------------------------------------
amn = readr::read_csv("data/Amniote_Database.csv") |>
  tidyr::unite(`Organism Name`, genus, species, sep = " ") |>
  dplyr::select(class, order, `Organism Name`, adult_body_mass_g, longevity_y, maximum_longevity_y) |>
  dplyr::filter(adult_body_mass_g >= 0 & longevity_y >= 0) |>
  dplyr::mutate(
    class = if_else(order == "Chiroptera", order, class),
    class = forcats::fct_relevel(class, "Mammalia", "Reptilia", "Aves", "Chiroptera"),
    volancy = dplyr::if_else(class %in% c("Aves", "Chiroptera"), class, "others")
  ) |>
  dplyr::arrange(class)

fit = lm(log10(longevity_y) ~ log10(adult_body_mass_g), data = amn)
# summary(fit)
amn_add_fit = modelr::add_predictions(amn, fit, var = "Pred") |>
  mutate(Residual = longevity_y - 10**Pred)


# Species selection ------------------------------------------------------------
df_ncbi = readr::read_tsv("data/amniota_ncbi_dataset.tsv") |>
  dplyr::rename(
    accession = `Assembly Accession`,
    genome_size = `Assembly Stats Total Sequence Length`,
    level = `Assembly Level`,
    N50 = `Assembly Stats Scaffold N50`,
    gene_count = `Annotation Count Gene Total`,
    protein_coding = `Annotation Count Gene Protein-coding`
  ) |>
  dplyr::filter(str_starts(accession, "GCF_")) |>
  #dplyr::filter(!accession == "GCF_004027225.2") |>  # remove Kakapo
  dplyr::select(
    accession, `Organism Name`, genome_size, level,
    N50, gene_count, protein_coding
  ) |>
  dplyr::inner_join(amn_add_fit, by = "Organism Name") |>
  dplyr::filter(level == "Chromosome" | order %in% c("Chiroptera", "Crocodilia"))

readr::write_tsv(df_ncbi, file = "out/species.tsv")
readr::write_lines(df_ncbi$`Organism Name`, "out/species.list")  # for TimeTree
readr::write_lines(df_ncbi$accession, "out/accession.list")
