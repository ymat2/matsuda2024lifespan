library(conflicted)
library(tidyverse)
library(RERconverge)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggtree)
library(ggimage)


mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

rer_results_all = readr::read_tsv("out/binary.nogap.tsv") |>
  #dplyr::left_join(aa_length, by = "symbol") |>
  #dplyr::filter(length >= 30) |>
  dplyr::left_join(sym2geneid, by = "symbol") |>
  tidyr::drop_na(gene_id)
rer_results = rer_results_all |>
  dplyr::filter(p.adj < 0.05) |>
  tibble::column_to_rownames(var = "gene_id") |>
  RERconverge::getStat() |>
  sort(decreasing=TRUE)


## Fig1a. Schematic tree

set.seed(43093)
bird.anc = 231
bat.anc = 175
species_tree = ape::read.tree("data/species.nwk") |>
  dplyr::as_tibble() |>
  ggtree::groupClade(c(bird.anc, bat.anc))
phylopic_info = data.frame(
  #node = c(bird.anc, bat.anc),
  node = c(102, 45),
  phylopic = c("c97400f5-fab6-4452-ab32-1ee071848f31",
               "18bfd2fc-f184-4c3a-b511-796aafcc70f6"),
  species = c("birds", "bats")
)
phylopic_info_acc = phylopic_info |> dplyr::mutate(node = c(68, 36))

#### accelerated genes

tr = species_tree |>
  dplyr::mutate(branch.length = dplyr::if_else(
    group == 1 & node != bird.anc,
    branch.length*5,
    branch.length*1)
  ) |>
  dplyr::mutate(branch.length = dplyr::if_else(
    group == 2 & node != bat.anc,
    branch.length*7,
    branch.length)
  ) |>
  ape::as.phylo() |>
  ggtree::groupClade(c(bird.anc, bat.anc))

pg = ggtree(tr) +
  aes(color = group) +
  scale_color_manual(values = c("#444444", "#D73027", "#D73027", "#D73027", "#D73027")) +
  theme(legend.position = "none")
tr_acc = pg %<+% phylopic_info_acc +
  geom_tiplab(aes(image = phylopic, color = species), geom = "phylopic", size = .1, offset = 200) +
  xlim(0, 700)

#### conserved genes

tr = species_tree |>
  dplyr::mutate(branch.length = dplyr::if_else(
    group == 1 | group == 2,
    branch.length*1/3,
    branch.length*1)
  ) |>
  ape::as.phylo() |>
  ggtree::groupClade(c(bird.anc, bat.anc))

pg = ggtree(tr) +
  aes(color = group) +
  scale_color_manual(values = c("#444444", "#4575B4", "#4575B4")) +
  theme(legend.position = "none")
tr_con = pg %<+% phylopic_info +
  geom_tiplab(aes(image = phylopic), geom = "phylopic", size = .1, offset = 200) +
  xlim(0, 700)

tr = cowplot::plot_grid(tr_acc, tr_con, ncol = 2, scale = .95)


## Fig.1b: Volcano plot

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
paml_accelerated_genes = paml_birds |>
  dplyr::inner_join(paml_bats, by = "symbol") |>
  dplyr::filter(p_adj_birds < 0.05 & p_adj_bats < 0.05 & diff_omega_birds > 0 & diff_omega_bats > 0)
paml_conserved_genes = paml_birds |>
  dplyr::inner_join(paml_bats, by = "symbol") |>
  dplyr::filter(p_adj_birds < 0.05 & p_adj_bats < 0.05 & diff_omega_birds < 0 & diff_omega_bats < 0)

df_volc = rer_results_all |>
  dplyr::mutate(type = dplyr::case_when(
    p.adj < 0.05 & Rho > 0 ~ "accelerated",
    p.adj < 0.05 & Rho < 0 ~ "conserved",
    .default = "not significant"
  )) |>
  dplyr::mutate(paml = dplyr::case_when(
    type == "accelerated" & symbol %in% paml_accelerated_genes$symbol ~ "paml_acc",
    type == "conserved" & symbol %in% paml_conserved_genes$symbol ~ "paml_cons",
    .default = "paml_other"
  )) |>
  dplyr::arrange(desc(paml))

vol = ggplot(df_volc) +
  aes(x = Rho, y = -log10(p.adj), color = paml, shape = paml, size = paml) +
  geom_point() +
  scale_color_manual(values = c("#D73027", "#2166ac", "#999999")) +
  scale_shape_manual(values = c(16, 16, 1)) +
  scale_size_manual(values = c(4, 4, 3)) +
  labs(x = "Rho", y = expression(paste("−", {log[10]}, "(adjusted p-value)", sep=""))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")


## Fig.1c: Dotplot

df_gse = readr::read_tsv("out/kegg_gsea_result.tsv") |>
  dplyr::select(Description, setSize, enrichmentScore, qvalue) |>
  dplyr::mutate(
    scoreType = -sign(enrichmentScore),
    enrichmentScore = abs(enrichmentScore),
    Description = forcats::fct_reorder(Description, enrichmentScore)
  )

dplot = ggplot(df_gse) +
  aes(x = Description, y = enrichmentScore, color = qvalue, size = setSize) +
  geom_point() +
  coord_flip() +
  scale_color_viridis_c(option = "E") +
  facet_grid(scoreType ~ ., scales = "free", space = "free") +
  labs(y = "Enrichment score", x = "", color = "q-value", size = "count") +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )


## Cowplot and save

fig_ab = cowplot::plot_grid(tr, vol, ncol = 2, labels = c("a", "b"), label_size = 36)
fig_abc = cowplot::plot_grid(fig_ab, dplot, nrow = 2)
ggsave(file = "images/evolutionary_rate_analysis.png", fig_abc, w = 12, h = 12, bg = "#ffffff")
