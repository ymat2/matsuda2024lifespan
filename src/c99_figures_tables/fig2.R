library(conflicted)
library(tidyverse)
library(RERconverge)
library(org.Hs.eg.db)


mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

kegg = clusterProfiler::download_KEGG("hsa")
lrp = kegg$KEGGPATHID2EXTID |>
  dplyr::filter(from == "hsa04211") |>
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

df4boxplot = dplyr::full_join(
  res_birds |>
    dplyr::filter(symbol %in% lrp$symbol) |>
    dplyr::rename(Rho_b = Rho, qval_b = p.adj) |>
    dplyr::select(symbol, Rho_b, qval_b),
  res_bats |>
    dplyr::filter(symbol %in% lrp$symbol) |>
    dplyr::rename(Rho_a = Rho, qval_a = p.adj) |>
    dplyr::select(symbol, Rho_a, qval_a),
  by = "symbol") |>
  dplyr::select(symbol, Rho_b, Rho_a) |>
  dplyr::rename("Birds" = Rho_b, "Bats" = Rho_a) |>
  dplyr::mutate(symbol = forcats::fct_inorder(symbol)) |>
  tidyr::pivot_longer(2:3, names_to = "group", values_to = "Rho") |>
  dplyr::left_join(lrp_sig, by = c("symbol", "group"))

birds_rho = df4boxplot |> dplyr::filter(group == "Birds")
bats_rho = df4boxplot |> dplyr::filter(group == "Bats")
wilcox = wilcox.test(birds_rho$Rho, bats_rho$Rho, paired = TRUE)

birds_rho_sig = df4boxplot |> dplyr::filter(group == "Birds" & sig == "Y")
bats_rho_sig = df4boxplot |> dplyr::filter(group == "Bats" & sig == "Y")
wilcox_sig = wilcox.test(birds_rho_sig$Rho, bats_rho_sig$Rho)

box = ggplot(df4boxplot) +
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
  labs(x = "", y = "Rho") +
  theme(legend.position = "none", panel.grid = element_blank())


df4tile = dplyr::full_join(
  res_birds |>
    dplyr::filter(symbol %in% lrp$symbol) |>
    dplyr::rename(Rho_b = Rho, qval_b = p.adj) |>
    dplyr::select(symbol, Rho_b, qval_b),
  res_bats |>
    dplyr::filter(symbol %in% lrp$symbol) |>
    dplyr::rename(Rho_a = Rho, qval_a = p.adj) |>
    dplyr::select(symbol, Rho_a, qval_a),
  by = "symbol") |>
  dplyr::filter(qval_a < 0.05 | qval_b < 0.05) |>
  dplyr::select(symbol, Rho_b, Rho_a) |>
  dplyr::rename("Birds" = Rho_b, "Bats" = "Rho_a") |>
  tidyr::pivot_longer(2:3, names_to = "group", values_to = "Rho")

tile = ggplot(df4tile) +
  aes(x = symbol, y = group, fill = Rho) +
  geom_tile() +
  coord_flip() +
  labs(x = "", y = "") +
  scale_fill_gradient(high = "#D73027", low = "#4575B4", na.value = "#999999") +
  theme_minimal(base_size = 20) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7, face = "bold"))


## oxygen

bird = rphylopic::get_phylopic(uuid = "eeb0df66-27ac-4cd9-852d-1af5a1f7679c")
bat = rphylopic::get_phylopic(uuid = "18bfd2fc-f184-4c3a-b511-796aafcc70f6")

ox_level = data.frame(
  Mya = c(252, 201, 190, 145, 66, 26, 0),
  o2 = c(20, 14, 11, 15, 19, 20, 21)
) |>
  ggplot() +
  geom_area(aes(x = Mya, y = o2), fill = "#DDDDDD") +
  scale_x_reverse(limits = c(252, -30)) +
  scale_y_continuous(limits = c(0, 35)) +
  annotate("text", x = 222.5, y = 10, size = 4, label = "Evolution of \nhypoxia tolerance",
           fontface = "bold", color = "#444444") +
  annotate("text", x = 122.5, y = 10, size = 4, label = "Adaptation to \nhigh oxygen level",
           fontface = "bold", color = "#444444") +
  geom_segment(
    aes(x = 250, y = 4, xend = 195, yend = 4),
    arrow = arrow(ends = "both", type = "open", length = unit(0.1, "inches")),
    color = "#444444") +
  geom_segment(
    aes(x = 70, y = 4, xend = 185, yend = 4),
    arrow = arrow(ends = "both", type = "open", length = unit(0.1, "inches")),
    color = "#444444") +
  geom_point(aes(150, 30), size = 5, color = "#009E73") +
  geom_segment(aes(x = 150, y = 30, xend = 0, yend = 30), color = "#009E73", linewidth = 2) +
  annotate("text", x = 200, y = 30, size = 6, label = "150 MYA", color = "#444444") +
  geom_point(aes(55, 25), size = 5, color = "#E69F00") +
  geom_segment(aes(x = 55, y = 25, xend = 0, yend = 25), color = "#E69F00", linewidth = 2) +
  annotate("text", x = 100, y = 25, size = 6, label = "50-60 MYA", color = "#444444") +
  rphylopic::add_phylopic(img = bird, x = -20, y = 32, ysize = 7, fill = "#999999") +
  rphylopic::add_phylopic(img = bat, x = -20, y = 25, ysize = 5, fill = "#999999") +
  labs(
    x = "Divergence time (MYA)",
    y = expression(paste({O[2]}, " level (%)"))
  ) +
  theme_classic(base_size = 20)


p_ab = cowplot::plot_grid(box, tile, labels = c("a", "b"), rel_widths = c(1, 1.33), label_size = 20)
p = cowplot::plot_grid(p_ab, ox_level, nrow=  2, labels = c("", "c"), label_size = 20)
ggsave(file = "images/longevity_regulating_pathway.png", p, w = 7, h = 7, bg = "#ffffff")
