library(conflicted)
library(tidyverse)


## Read data

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

lht = readr::read_tsv("out/species.tsv") |>
  dplyr::mutate(
    rownm = stringr::str_replace(`Organism Name`, " ", "_"),
    class = factor(class, levels = c("Chiroptera", "Aves", "Mammalia", "Reptilia")),
  ) |>
  tibble::column_to_rownames(var = "rownm")

## Statistical analysis

residual_birds_all = amn_add_fit |> dplyr::filter(volancy == "Aves")
residual_bats_all = amn_add_fit |> dplyr::filter(volancy == "Chiroptera")
residual_others_all = amn_add_fit |> dplyr::filter(volancy == "others")

wilcox_birds_all = wilcox.test(
  x = residual_birds_all$Residual,
  y = residual_others_all$Residual,
  paired = FALSE
)

wilcox_bats_all = wilcox.test(
  x = residual_bats_all$Residual,
  y = residual_others_all$Residual,
  paired = FALSE
)


residual_birds_use = lht |> dplyr::filter(volancy == "Aves")
residual_bats_use = lht |> dplyr::filter(volancy == "Chiroptera")
residual_others_use = lht |> dplyr::filter(volancy == "others")

wilcox_birds_use = wilcox.test(
  x = residual_birds_use$Residual,
  y = residual_others_use$Residual,
  paired = FALSE
)

wilcox_bats_use = wilcox.test(
  x = residual_bats_use$Residual,
  y = residual_others_use$Residual,
  paired = FALSE
)


## Plot

plot_all = ggplot(amn_add_fit) +
  aes(x = adult_body_mass_g, y = longevity_y) +
  geom_point(aes(color = volancy, shape = volancy), size = 2) +
  geom_line(aes(y = 10**Pred), color = "#444444", linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("#009E73", "#E69F00", "#CCCCCC")) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10() +
  labs(x = "Body mass (g)", y = "Longevity (years)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 0),
    legend.background = element_rect(color = "#444444", linewidth = 1),
    legend.title = element_blank(),
    panel.grid = element_blank()
    )

box_all = ggplot(amn_add_fit) +
  aes(x = volancy, y = Residual, color = volancy) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = volancy), height = 0, width = .25, alpha = .5) +
  geom_segment(x = 1, xend = 3, y = 90, yend = 90, color = "#444444") +
  annotate("text", x = 2, y = 100, label = "< 1.0e-10", size = 5, color = "#444444") +
  geom_segment(x = 2, xend = 3, y = 60, yend = 60, color = "#444444") +
  annotate("text", x = 2.5, y = 70, label = "< 1.0e-10", size = 5, color = "#444444") +
  scale_color_manual(values = c("#009E73", "#E69F00", "#999999")) +
  labs(x = "", y = "Residual longevity (years)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
    )


plot_use = ggplot(lht) +
  aes(x = adult_body_mass_g, y = longevity_y) +
  geom_point(aes(color = volancy, shape = volancy), size = 3) +
  geom_line(aes(y = 10**Pred), color = "#444444", linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("#009E73", "#E69F00", "#CCCCCC")) +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10() +
  labs(x = "Body mass (g)", y = "Longevity (years)", color = "", shape = "") +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", panel.grid = element_blank())

box_use = ggplot(lht) +
  aes(x = volancy, y = Residual, color = volancy) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = volancy), height = 0, width = .25, alpha = .5) +
  geom_segment(x = 1, xend = 3, y = 90, yend = 90, color = "#444444") +
  annotate("text", x = 2, y = 100, size = 5, color = "#444444",
           label = wilcox_birds_use$p.value |> round(digits = 3)) +
  geom_segment(x = 2, xend = 3, y = 60, yend = 60, color = "#444444") +
  annotate("text", x = 2.5, y = 70, size = 5, color = "#444444",
           label = wilcox_bats_use$p.value |> round(digits = 4)) +
  scale_color_manual(values = c("#009E73", "#E69F00", "#999999")) +
  labs(x = "", y = "Residual longevity (years)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

figs1 = cowplot::plot_grid(plot_all, box_all, plot_use, box_use, nrow = 2, scale = .95,
                   rel_widths = c(1.41, 1), align = "hv", axis = "tb",
                   labels = c("a", "", "b", ""), label_size = 22)
ggsave("images/lifehistory_trait_distribution.png", figs1, w = 9, h = 9, bg = "#FFFFFF")
