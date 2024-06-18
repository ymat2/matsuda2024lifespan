library(conflicted)
library(tidyverse)
library(org.Hs.eg.db)
library(ggvenn)


### Reading data

mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

rer_results = readr::read_tsv("out/binary.nogap.tsv") |>
  #dplyr::left_join(aa_length, by = "symbol") |>
  #dplyr::filter(length >= 30) |>
  dplyr::left_join(sym2geneid, by = "symbol") |>
  tidyr::drop_na(gene_id)

paml_birds = readr::read_tsv("out/paml_aves_result.tsv") |>
  dplyr::mutate(
    symbol = stringr::str_split(symbol, "/", simplify = TRUE)[, 2],
    p_adj_birds = stats::p.adjust(`p-value`, method = "BH"),
    diff_omega_birds = `foreground-omega` - `background-omega`
  ) |>
  dplyr::select(symbol, p_adj_birds, diff_omega_birds)
paml_bats = readr::read_tsv("out/paml_chiroptera_result.tsv") |>
  dplyr::mutate(
    symbol = stringr::str_split(symbol, "/", simplify = TRUE)[, 2],
    p_adj_bats = stats::p.adjust(`p-value`, method = "BH"),
    diff_omega_bats = `foreground-omega` - `background-omega`
  ) |>
  dplyr::select(symbol, p_adj_bats, diff_omega_bats)

paml_results = dplyr::inner_join(paml_birds, paml_bats, by = "symbol")


### Permutation test

permt = function(N1, pop1, N2, pop2, overlap, perm = 1000) {
  n = 0
  L = numeric(perm)
  for (i in 1:perm) {
    subpop1 = base::sample(pop1, N1)
    subpop2 = base::sample(pop2, N2)
    l = length(base::intersect(subpop1, subpop2))
    L[i] = l
    if (l >= overlap) n = n + 1
  }

  df4plot = base::data.frame("L" = L)
  p = ggplot(df4plot) +
    aes(x = L) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = overlap, color = "#D73027") +
    labs(x = "Number of overlapping genes", y = "Frequency") +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank())

  plot(p)

  pval = n/perm

  return(list(pvalue = pval, plot = p))
}

### Accelerated genes

rer_accelerated_genes = rer_results |>
  dplyr::filter(p.adj < 0.05 & Rho > 0)  # 2042

paml_accelerated_genes = paml_results |>
  dplyr::filter(p_adj_birds < 0.05 & diff_omega_birds > 0) |>
  dplyr::filter(p_adj_bats < 0.05 & diff_omega_bats > 0)

accelerated_gene_list = list(
  RERconverge = rer_accelerated_genes$symbol,
  PAML = paml_accelerated_genes$symbol
)

va = ggvenn::ggvenn(
  accelerated_gene_list,
  fill_color = c("#FFFFFF", "#FFFFFF"),
  fill_alpha = 1,
  stroke_color = "#444444",
  set_name_color = "#444444",
  set_name_size = 6,
  text_color = "#444444",
  text_size = 4) +
  scale_y_continuous(expand = c(.1, .1))
test_acc = permt(length(rer_accelerated_genes$symbol), pop1 = rer_results$symbol,
                 length(paml_accelerated_genes$symbol), pop2 = paml_results$symbol,
                 overlap = 397, perm = 5000)
ba = test_acc$plot


### Conserved genes

rer_conserved_genes = rer_results |>
  dplyr::filter(p.adj < 0.05 & Rho < 0)  # 1337

paml_conserved_genes = dplyr::inner_join(paml_birds, paml_bats, by = "symbol") |>
  dplyr::filter(p_adj_birds < 0.05 & diff_omega_birds < 0) |>
  dplyr::filter(p_adj_bats < 0.05 & diff_omega_bats < 0)

conserved_gene_list = list(
  RERconverge = rer_conserved_genes$symbol,
  PAML = paml_conserved_genes$symbol
)

vc = ggvenn::ggvenn(
  conserved_gene_list,
  fill_color = c("#FFFFFF", "#FFFFFF"),
  fill_alpha = 1,
  stroke_color = "#444444",
  set_name_color = "#444444",
  set_name_size = 6,
  text_color = "#444444",
  text_size = 4) +
  scale_y_continuous(expand = c(.1, .1))
test_cons = permt(length(paml_conserved_genes$symbol), pop1 = rer_results$symbol,
                  length(paml_conserved_genes$symbol), pop2 = paml_results$symbol,
                  overlap = 116, perm = 5000)
bc = test_cons$plot


### cowplot

p = cowplot::plot_grid(va, ba, vc, bc, nrow = 2, scale = .99, labels = letters[1:4])
ggsave(file = "images/overlap_with_permtation.png", p, w = 8, h = 8, bg = "#ffffff")
