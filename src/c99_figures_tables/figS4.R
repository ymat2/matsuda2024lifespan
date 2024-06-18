library(conflicted)
library(tidyverse)
library(org.Hs.eg.db)


cancer_genes = readr::read_tsv("data/cancerGeneList.tsv") |>
  dplyr::rename(symbol = `Hugo Symbol`) |>
  dplyr::select(symbol, `Is Oncogene`, `Is Tumor Suppressor Gene`)
oncogenes = cancer_genes |>
  dplyr::filter(`Is Oncogene` == "Yes")
tsgenes = cancer_genes |>
  dplyr::filter(`Is Tumor Suppressor Gene` == "Yes")

mapped_genes = mappedkeys(org.Hs.egSYMBOL2EG)
sym2geneid = as.data.frame(org.Hs.egSYMBOL2EG[mapped_genes])

rer_results = readr::read_tsv("out/binary.nogap.tsv") |>
  dplyr::left_join(sym2geneid, by = "symbol") |>
  tidyr::drop_na(gene_id)


oncogenes = oncogenes |>
  dplyr::inner_join(rer_results, by = "symbol") |>
  dplyr::mutate(type = "oncogene") |>
  dplyr::select(symbol, Rho, p.adj, type)

tsgenes = tsgenes |>
  dplyr::inner_join(rer_results, by = "symbol") |>
  dplyr::mutate(type = "Tumor suppressor gene") |>
  dplyr::select(symbol, Rho, p.adj, type)

wilcox_onco = wilcox.test(oncogenes$Rho, rer_results$Rho, paired = FALSE)
wilcox_tsg = wilcox.test(tsgenes$Rho, rer_results$Rho, paired = FALSE)

df4boxplot = dplyr::bind_rows(oncogenes, tsgenes)
p = ggplot(df4boxplot) +
  aes(x = type, y = Rho) +
  geom_boxplot(outliers = FALSE, linewidth = 1) +
  geom_jitter(height = 0, width = .2, size = 3, alpha = .5, shape = 16) +
  geom_hline(aes(yintercept = 0), color = "darkred", linetype = "dashed", linewidth = 1.5) +
  annotate("text", x = 1, y = .6, label = "P = 0.004", size = 5, color = "#444444") +
  annotate("text", x = 2, y = .6, label = "P = 0.097", size = 5, color = "#444444") +
  scale_y_continuous(limits = c(NA, 0.65)) +
  theme_bw(base_size = 20) +
  #labs(x = "", y = "", title = "Longevity regulating pathway\n- multiple species") +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(file = "images/rho_cancer_genes.png", p, w = 7, h = 4)
