library(conflicted)
library(tidyverse)
library(RERconverge)
library(org.Hs.eg.db)


genage_genes = readr::read_csv("data/genage_human.csv")

mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

kegg = clusterProfiler::download_KEGG("hsa")
lrp_mlt = kegg$KEGGPATHID2EXTID |>
  dplyr::filter(from == "hsa04213") |>  # multiple species
  dplyr::rename(gene_id = to) |>
  dplyr::left_join(sym2geneid, by = "gene_id")

trees_nogap = readr::read_rds("data/trees.nogap.rds")
RERs = readr::read_rds("data/RER_all.nogap.rds")

# Binary trait analysis --------------------------------------------------------

lht = readr::read_tsv("out/species.tsv") |>
  dplyr::mutate(rownm = stringr::str_replace(`Organism Name`, " ", "_")) |>
  tibble::column_to_rownames(var = "rownm")


### Aves -----------------------------------------------------------------------

tree = foreground2Tree(
  rownames(dplyr::filter(lht, class == "Aves")),
  trees_nogap,
  clade="all",  # alternatively [all, terminal, ancestral]
  useSpecies=rownames(lht)
)
phenv=tree2Paths(tree, trees_nogap)

res_birds = RERconverge::correlateWithBinaryPhenotype(
  RERs,
  phenv,
  min.sp=10,
  min.pos=10,
  weighted="auto"
) |>
  tibble::rownames_to_column("symbol") |>
  dplyr::mutate(symbol = str_replace(symbol, ".pep.aln.*", "")) |>
  dplyr::arrange(P) |>
  tidyr::drop_na()
readr::write_tsv(res_birds, file = "out/binary.birds.tsv")

### Chiroptera -----------------------------------------------------------------

tree = foreground2Tree(
  rownames(dplyr::filter(lht, class == "Chiroptera")),
  trees_nogap,
  clade="all",  # alternatively [all, terminal, ancestral]
  useSpecies=rownames(lht)
)
phenv=tree2Paths(tree, trees_nogap)

res_bats = RERconverge::correlateWithBinaryPhenotype(
  RERs,
  phenv,
  min.sp=10,
  min.pos=10,
  weighted="auto"
) |>
  tibble::rownames_to_column("symbol") |>
  dplyr::mutate(symbol = str_replace(symbol, ".pep.aln.*", "")) |>
  dplyr::arrange(P) |>
  tidyr::drop_na()
readr::write_tsv(res_bats, file = "out/binary.bats.tsv")

lrp_sig = dplyr::bind_rows(
  res_birds |> dplyr::mutate(group = "Birds"),
  res_bats |> dplyr::mutate(group = "Bats")
) |>
  dplyr::mutate(sig = dplyr:::if_else(p.adj < 0.05, "Y", "N")) |>
  dplyr::select(symbol, group, sig)

## Boxplot for all genes -------------------------------------------------------

df4boxplot_all = dplyr::full_join(
  res_birds |>
    dplyr::rename(Rho_b = Rho, qval_b = p.adj) |>
    dplyr::select(symbol, Rho_b, qval_b),
  res_bats |>
    dplyr::rename(Rho_a = Rho, qval_a = p.adj) |>
    dplyr::select(symbol, Rho_a, qval_a),
  by = "symbol") |>
  dplyr::select(symbol, Rho_b, Rho_a) |>
  dplyr::rename("Birds" = Rho_b, "Bats" = Rho_a) |>
  dplyr::mutate(symbol = forcats::fct_inorder(symbol)) |>
  tidyr::pivot_longer(2:3, names_to = "group", values_to = "Rho")

birds_rho = df4boxplot_all |> dplyr::filter(group == "Birds")
bats_rho = df4boxplot_all |> dplyr::filter(group == "Bats")
wilcox = wilcox.test(birds_rho$Rho, bats_rho$Rho, alternative = "less", paired = TRUE)

box_all = ggplot(df4boxplot_all) +
  aes(x = group, y = Rho, color = group) +
  geom_hline(aes(yintercept = 0), color = "#444444", linetype = "solid", linewidth = 1.5) +
  geom_jitter(height = 0, width = .2, size = 3, alpha = 1) +
  geom_boxplot(outliers = FALSE, linewidth = 1, shape = 16) +
  geom_segment(x = 1, xend = 2, y = .6, yend = .6, color = "#444444") +
  annotate("text", x = 1.5, y = .65, label = "N.S.", size = 5, color = "#444444") +
  scale_color_manual(values = c("#f5d899", "#99d8c7")) +
  scale_y_continuous(limits = c(NA, 0.7)) +
  theme_bw(base_size = 20) +
  labs(x = "", y = "Rho", title = "All genes") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(color = "#444444", hjust = .5, size = 14)
  )


## Boxplot for lrp multiple species --------------------------------------------

df4boxplot_mlt = dplyr::full_join(
  res_birds |>
    dplyr::filter(symbol %in% lrp_mlt$symbol) |>
    dplyr::rename(Rho_b = Rho, qval_b = p.adj) |>
    dplyr::select(symbol, Rho_b, qval_b),
  res_bats |>
    dplyr::filter(symbol %in% lrp_mlt$symbol) |>
    dplyr::rename(Rho_a = Rho, qval_a = p.adj) |>
    dplyr::select(symbol, Rho_a, qval_a),
  by = "symbol") |>
  dplyr::select(symbol, Rho_b, Rho_a) |>
  dplyr::rename("Birds" = Rho_b, "Bats" = Rho_a) |>
  dplyr::mutate(symbol = forcats::fct_inorder(symbol)) |>
  tidyr::pivot_longer(2:3, names_to = "group", values_to = "Rho") |>
  dplyr::left_join(lrp_sig, by = c("symbol", "group"))

birds_rho = df4boxplot_mlt |> dplyr::filter(group == "Birds")
bats_rho = df4boxplot_mlt |> dplyr::filter(group == "Bats")
wilcox = wilcox.test(birds_rho$Rho, bats_rho$Rho, paired = TRUE)

box_mlt = ggplot(df4boxplot_mlt) +
  aes(x = group, y = Rho, color = group) +
  geom_hline(aes(yintercept = 0), color = "#444444", linetype = "solid", linewidth = 1.5) +
  geom_boxplot(outliers = FALSE, linewidth = 1, shape = 16) +
  geom_jitter(aes(shape = sig), height = 0, width = .2, size = 3, alpha = .7) +
  geom_segment(x = 1, xend = 2, y = .45, yend = .45, color = "#444444") +
  annotate("text", x = 1.5, y = .5, label = "<0.001", size = 5, color = "#444444") +
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  scale_shape_manual(values = c("Y" = 16, "N" = 1)) +
  scale_y_continuous(limits = c(NA, 0.55)) +
  theme_bw(base_size = 20) +
  labs(x = "", y = "", title = "Longevity regulating pathway\n- multiple species") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(color = "#444444", hjust = .5, size = 14)
    )


## Boxplot for geneage ---------------------------------------------------------

df4boxplot_genage = dplyr::full_join(
  res_birds |>
    dplyr::filter(symbol %in% genage_genes$symbol) |>
    dplyr::rename(Rho_b = Rho, qval_b = p.adj) |>
    dplyr::select(symbol, Rho_b, qval_b),
  res_bats |>
    dplyr::filter(symbol %in% genage_genes$symbol) |>
    dplyr::rename(Rho_a = Rho, qval_a = p.adj) |>
    dplyr::select(symbol, Rho_a, qval_a),
  by = "symbol") |>
  dplyr::select(symbol, Rho_b, Rho_a) |>
  dplyr::rename("Birds" = Rho_b, "Bats" = Rho_a) |>
  dplyr::mutate(symbol = forcats::fct_inorder(symbol)) |>
  tidyr::pivot_longer(2:3, names_to = "group", values_to = "Rho") |>
  dplyr::left_join(lrp_sig, by = c("symbol", "group"))

birds_rho = df4boxplot_genage |> dplyr::filter(group == "Birds")
bats_rho = df4boxplot_genage |> dplyr::filter(group == "Bats")
wilcox = wilcox.test(birds_rho$Rho, bats_rho$Rho, paired = TRUE)

box_genage = ggplot(df4boxplot_genage) +
  aes(x = group, y = Rho, color = group) +
  geom_hline(aes(yintercept = 0), color = "#444444", linetype = "solid", linewidth = 1.5) +
  geom_boxplot(outliers = FALSE, linewidth = 1, shape = 16) +
  geom_jitter(aes(shape = sig), height = 0, width = .2, size = 3, alpha = .7) +
  geom_segment(x = 1, xend = 2, y = .45, yend = .45, color = "#444444") +
  annotate("text", x = 1.5, y = .5, label = "N.S.", size = 5, color = "#444444") +
  scale_color_manual(values = c("#E69F00", "#009E73")) +
  scale_shape_manual(values = c("Y" = 16, "N" = 1)) +
  scale_y_continuous(limits = c(NA, 0.55)) +
  theme_bw(base_size = 20) +
  labs(x = "", y = "", title = "GenAge") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(color = "#444444", hjust = .5, size = 14)
  )


## Cowplot ---------------------------------------------------------------------

p = cowplot::plot_grid(box_all, box_mlt, box_genage, nrow = 1, labels = c("a", "b", "c"),
                       label_size = 18, align = "hv", axis = "tb")
ggsave(file = "images/boxplot_other_genes.png", p, w = 11, h = 5, bg = "#FFFFFF")
