library(conflicted)
library(tidyverse)


gene2desc = readr::read_tsv("data/symbol2description.tsv") |>
  dplyr::rename(symbol = `Gene name`, description = `Gene description`) |>
  dplyr::mutate(description = stringr::str_remove(description, " \\[.*\\]"))


csubst_result = readr::read_tsv("data/csubst_result.tsv", show_col_types = FALSE) |>
  dplyr::select(symbol, omegaCany2spe, OCNany2spe) |>
  dplyr::left_join(gene2desc, by = "symbol") |>
  dplyr::relocate(description, .after = symbol) |>
  dplyr::filter(omegaCany2spe > 3 & OCNany2spe > 3) |>
  dplyr::arrange(desc(omegaCany2spe)) |>
  dplyr::rename("Gene" = symbol, "Description" = description,  "$\\omega_C$" = omegaCany2spe, "${O_C}^N$" = OCNany2spe)

knitr::kable(csubst_result, align = "llrr")
