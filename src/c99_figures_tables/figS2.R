library(conflicted)
library(tidyverse)
library(org.Hs.eg.db)


mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

rer_results = readr::read_tsv("out/binary.nogap.tsv") |>
  #dplyr::left_join(aa_length, by = "symbol") |>
  #dplyr::filter(length >= 30) |>
  dplyr::left_join(sym2geneid, by = "symbol") |>
  #dplyr::filter(N > 200) |>
  tidyr::drop_na(gene_id)

pjit = ggplot(rer_results) +
  aes(x = N, y = -log10(p.adj)) +
  geom_point(size = 3, alpha = .6, shape = 16) +
  #geom_smooth(method = "lm", formula = y~x, color = "darkred") +
  labs(
    x = "Number of branches",
    y = expression(paste({-log[10]}, "(adjusted P)", sep=""))
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())


nbranch = rer_results$N |> unique() |> sort()
psig = numeric(length(nbranch))
for (i in 1:length(nbranch)) {
  df = rer_results |>
    dplyr::filter(N >= nbranch[i]) |>
    dplyr::mutate(qval = p.adjust(P, method = "BH"))
  dfsig = df |> dplyr::filter(qval < 0.05)
  nsig = nrow(dfsig)
  psig[i] = nsig
}
pcol = ggplot(data.frame("N" = nbranch, "Psig" = psig)) +
  aes(x = N, y = Psig) +
  geom_line() +
  labs(
    x = "Minimum number of branches of genes",
    y = "Number of significant genes"
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

p = cowplot::plot_grid(pjit, pcol, nrow = 1, labels = c("a", "b"))
ggsave(file = "images/number_of_branch.png", p, w = 10, h = 4)
