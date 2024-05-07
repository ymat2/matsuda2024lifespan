library(conflicted)
library(tidyverse)
library(RERconverge)
library(org.Hs.eg.db)
library(clusterProfiler)


mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])


rer_results_all = readr::read_tsv("out/binary.nogap.tsv") |>
  dplyr::left_join(sym2geneid, by = "symbol") |>
  tidyr::drop_na(gene_id)
rer_results = rer_results_all |>
  dplyr::filter(p.adj < 0.05) |>
  tibble::column_to_rownames(var = "gene_id") |>
  RERconverge::getStat() |>
  sort(decreasing=TRUE)

## GSE for accelerated/conserved genes ----------------------------------------

kse = clusterProfiler::gseKEGG(
  geneList = rer_results,
  organism = "hsa",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = TRUE) |>
  clusterProfiler::setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
kegg_result = kse@result |> dplyr::filter(qvalue < 0.05)
readr::write_tsv(kegg_result, file = "out/kegg_gsea_result.tsv")
