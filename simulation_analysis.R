library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)
library(grid)

# read file ---- ###############################################################

sim_data_all <- fread(file = "prevCI_hypgeom_simData.csv")

# N50 ---- ####################################################################

sim_data_N50 <- sim_data_all %>% filter(N == 50)

df_cov_N50 <- sim_data_N50 %>%
  pivot_longer(
    cols = c(coverage_exact_B_mean, coverage_simul_B_mean, coverage_Bayes_mean,
             coverage_FPC_mean, coverage_Ge_mean, coverage_binom_B_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

method_colors <- c("exact_B" = "#E41A1C",
                   "simul_B" = "#4CBB17",
                   "Bayes" = "#377EB8",
                   "FPC" = "#E41A1C",
                   "Ge" = "#377EB8",
                   "binom_B" = "#000000")

method_linetypes <- c("exact_B" = "1",
                      "simul_B" = "1",
                      "Bayes" = "1",
                      "FPC" = "2",
                      "Ge" = "2",
                      "binom_B" = "3")

df_cov_N50 <- df_cov_N50 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_cov_50_09 <-
  ggplot(df_cov_N50 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N50, Se == 0.9 & Sp == 0.9)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Mean coverage propability",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_50_099 <-
  ggplot(df_cov_N50 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N50, Se == 0.99 & Sp == 0.99)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N50 <- sim_data_N50 %>%
  pivot_longer(
    cols = c(width_exact_B_mean, width_simul_B_mean, width_Bayes_mean,
             width_FPC_mean, width_Ge_mean, width_binom_B_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

df_wid_N50 <- df_wid_N50 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_wid_50_09 <-
  ggplot(df_wid_N50 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N50, Se == 0.9 & Sp == 0.9)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Interval width",
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

plot_wid_50_099 <-
  ggplot(df_wid_N50 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N50, Se == 0.99 & Sp == 0.99)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

df_loc_N50 <- sim_data_N50 %>%
  pivot_longer(
    cols = c(loc_exact_mean, loc_simul_mean, loc_Bayes_mean,
             loc_Ge_mean, loc_RG_mean),
    names_to = "method",
    names_pattern = "loc_(.*)_mean",
    values_to = "loc_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact",
                                            "simul",
                                            "Bayes",
                                            "Ge",
                                            "RG"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, loc_mean))

method_colors_loc <- c("exact" = "#E41A1C",
                       "simul" = "#4CBB17",
                       "Bayes" = "#377EB8",
                       "Ge" = "#377EB8",
                       "RG" = "#000000")

method_linetypes_loc <- c("exact" = "1",
                          "simul" = "1",
                          "Bayes" = "1",
                          "Ge" = "2",
                          "RG" = "3")

df_loc_N50 <- df_loc_N50 %>%
  mutate(linetype = method_linetypes_loc[as.character(method)])

plot_loc_50_09 <-
  ggplot(df_loc_N50 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = "Point estimate",
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_loc_50_099 <-
  ggplot(df_loc_N50 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

legend_plot <-
  ggplot(df_wid_N50,
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Simul hypgeom CI",
                       "Bayes hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)",
                       "Binom CI (Reiczigel et al. 2010) / Rogan-Gladen estimate"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 2, 2, 2)), nrow = 1))

legend <- cowplot::get_legend(legend_plot)

title_grob <- textGrob("Population size: 50 (sample size in panels)\n",
                       gp = gpar(fontsize = 15, fontface = "bold"))

subheadings <- arrangeGrob(
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.9), gp = gpar(fontsize = 12)),
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.99), gp = gpar(fontsize = 12)),
  ncol = 2
)

plot_N50 <- 
grid.arrange(subheadings,
             arrangeGrob(plot_cov_50_09, plot_cov_50_099,
                         plot_wid_50_09, plot_wid_50_099,
                         plot_loc_50_09, plot_loc_50_099,
                         ncol = 2),
             legend,
             ncol = 1,
             heights = c(0.5, 10, 1),
             top = title_grob)

ggsave("plot_N50.svg", plot = plot_N50, width = 12.8, height = 10.6, units = "in", dpi = 600)
ggsave("plot_N50.pdf", plot = plot_N50, width = 12.8, height = 10.6, units = "in", dpi = 600)

# N100 ---- ####################################################################

sim_data_N100 <- sim_data_all %>% filter(N == 100)

df_cov_N100 <- sim_data_N100 %>%
  pivot_longer(
    cols = c(coverage_exact_B_mean, coverage_simul_B_mean, coverage_Bayes_mean,
             coverage_FPC_mean, coverage_Ge_mean, coverage_binom_B_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

method_colors <- c("exact_B" = "#E41A1C",
                   "simul_B" = "#4CBB17",
                   "Bayes" = "#377EB8",
                   "FPC" = "#E41A1C",
                   "Ge" = "#377EB8",
                   "binom_B" = "#000000")

method_linetypes <- c("exact_B" = "1",
                      "simul_B" = "1",
                      "Bayes" = "1",
                      "FPC" = "2",
                      "Ge" = "2",
                      "binom_B" = "3")

df_cov_N100 <- df_cov_N100 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_cov_100_09 <-
  ggplot(df_cov_N100 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N100, Se == 0.9 & Sp == 0.9)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Mean coverage propability",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_100_099 <-
  ggplot(df_cov_N100 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N100, Se == 0.99 & Sp == 0.99)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N100 <- sim_data_N100 %>%
  pivot_longer(
    cols = c(width_exact_B_mean, width_simul_B_mean, width_Bayes_mean,
             width_FPC_mean, width_Ge_mean, width_binom_B_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

df_wid_N100 <- df_wid_N100 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_wid_100_09 <-
  ggplot(df_wid_N100 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N100, Se == 0.9 & Sp == 0.9)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Interval width",
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

plot_wid_100_099 <-
  ggplot(df_wid_N100 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N100, Se == 0.99 & Sp == 0.99)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

df_loc_N100 <- sim_data_N100 %>%
  pivot_longer(
    cols = c(loc_exact_mean, loc_simul_mean, loc_Bayes_mean,
             loc_Ge_mean, loc_RG_mean),
    names_to = "method",
    names_pattern = "loc_(.*)_mean",
    values_to = "loc_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact",
                                            "simul",
                                            "Bayes",
                                            "Ge",
                                            "RG"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, loc_mean))

method_colors_loc <- c("exact" = "#E41A1C",
                       "simul" = "#4CBB17",
                       "Bayes" = "#377EB8",
                       "Ge" = "#377EB8",
                       "RG" = "#000000")

method_linetypes_loc <- c("exact" = "1",
                          "simul" = "1",
                          "Bayes" = "1",
                          "Ge" = "2",
                          "RG" = "3")

df_loc_N100 <- df_loc_N100 %>%
  mutate(linetype = method_linetypes_loc[as.character(method)])

plot_loc_100_09 <-
  ggplot(df_loc_N100 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = "Point estimate",
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_loc_100_099 <-
  ggplot(df_loc_N100 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

legend_plot <-
  ggplot(df_wid_N100,
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Simul hypgeom CI",
                       "Bayes hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)",
                       "Binom CI (Reiczigel et al. 2010) / Rogan-Gladen estimate"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 2, 2, 2)), nrow = 1))

legend <- cowplot::get_legend(legend_plot)

title_grob <- textGrob("Population size: 100 (sample size in panels)\n",
                       gp = gpar(fontsize = 15, fontface = "bold"))

subheadings <- arrangeGrob(
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.9), gp = gpar(fontsize = 12)),
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.99), gp = gpar(fontsize = 12)),
  ncol = 2
)

plot_N100 <- 
  grid.arrange(subheadings,
               arrangeGrob(plot_cov_100_09, plot_cov_100_099,
                           plot_wid_100_09, plot_wid_100_099,
                           plot_loc_100_09, plot_loc_100_099,
                           ncol = 2),
               legend,
               ncol = 1,
               heights = c(0.5, 10, 1),
               top = title_grob)

ggsave("plot_N100.svg", plot = plot_N100, width = 12.8, height = 10.6, units = "in", dpi = 600)
ggsave("plot_N100.pdf", plot = plot_N100, width = 12.8, height = 10.6, units = "in", dpi = 600)

# N200 ---- ####################################################################

sim_data_N200 <- sim_data_all %>% filter(N == 200)

df_cov_N200 <- sim_data_N200 %>%
  pivot_longer(
    cols = c(coverage_exact_B_mean, coverage_simul_B_mean, coverage_Bayes_mean,
             coverage_FPC_mean, coverage_Ge_mean, coverage_binom_B_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

method_colors <- c("exact_B" = "#E41A1C",
                   "simul_B" = "#4CBB17",
                   "Bayes" = "#377EB8",
                   "FPC" = "#E41A1C",
                   "Ge" = "#377EB8",
                   "binom_B" = "#000000")

method_linetypes <- c("exact_B" = "1",
                      "simul_B" = "1",
                      "Bayes" = "1",
                      "FPC" = "2",
                      "Ge" = "2",
                      "binom_B" = "3")

df_cov_N200 <- df_cov_N200 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_cov_200_09 <-
  ggplot(df_cov_N200 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N200, Se == 0.9 & Sp == 0.9)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Mean coverage propability",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_200_099 <-
  ggplot(df_cov_N200 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  coord_cartesian(ylim = c(max(0.9, min(subset(df_cov_N200, Se == 0.99 & Sp == 0.99)$coverage_mean)), 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N200 <- sim_data_N200 %>%
  pivot_longer(
    cols = c(width_exact_B_mean, width_simul_B_mean, width_Bayes_mean,
             width_FPC_mean, width_Ge_mean, width_binom_B_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

df_wid_N200 <- df_wid_N200 %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_wid_200_09 <-
  ggplot(df_wid_N200 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N200, Se == 0.9 & Sp == 0.9)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Interval width",
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

plot_wid_200_099 <-
  ggplot(df_wid_N200 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  geom_hline(yintercept = 1, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  coord_cartesian(ylim = c(0, min(1.5, max(subset(df_wid_N200, Se == 0.99 & Sp == 0.99)$width_mean)))) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = NULL,
       x = "True prevalence") +
  theme_bw() +
  theme(legend.position = "none")

df_loc_N200 <- sim_data_N200 %>%
  pivot_longer(
    cols = c(loc_exact_mean, loc_simul_mean, loc_Bayes_mean,
             loc_Ge_mean, loc_RG_mean),
    names_to = "method",
    names_pattern = "loc_(.*)_mean",
    values_to = "loc_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact",
                                            "simul",
                                            "Bayes",
                                            "Ge",
                                            "RG"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, loc_mean))

method_colors_loc <- c("exact" = "#E41A1C",
                       "simul" = "#4CBB17",
                       "Bayes" = "#377EB8",
                       "Ge" = "#377EB8",
                       "RG" = "#000000")

method_linetypes_loc <- c("exact" = "1",
                          "simul" = "1",
                          "Bayes" = "1",
                          "Ge" = "2",
                          "RG" = "3")

df_loc_N200 <- df_loc_N200 %>%
  mutate(linetype = method_linetypes_loc[as.character(method)])

plot_loc_200_09 <-
  ggplot(df_loc_N200 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = "Point estimate",
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_loc_200_099 <-
  ggplot(df_loc_N200 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = loc_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors_loc) +
  labs(y = NULL,
       x = "True prevalence") +
  geom_abline(slope = 1, intercept = 0, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

legend_plot <-
  ggplot(df_wid_N200,
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Simul hypgeom CI",
                       "Bayes hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)",
                       "Binom CI (Reiczigel et al. 2010) / Rogan-Gladen estimate"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 2, 2, 2)), nrow = 1))

legend <- cowplot::get_legend(legend_plot)

title_grob <- textGrob("Population size: 200 (sample size in panels)\n",
                       gp = gpar(fontsize = 15, fontface = "bold"))

subheadings <- arrangeGrob(
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.9), gp = gpar(fontsize = 12)),
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.99), gp = gpar(fontsize = 12)),
  ncol = 2
)

plot_N200 <- 
  grid.arrange(subheadings,
               arrangeGrob(plot_cov_200_09, plot_cov_200_099,
                           plot_wid_200_09, plot_wid_200_099,
                           plot_loc_200_09, plot_loc_200_099,
                           ncol = 2),
               legend,
               ncol = 1,
               heights = c(0.5, 10, 1),
               top = title_grob)

ggsave("plot_N200.svg", plot = plot_N200, width = 12.8, height = 10.6, units = "in", dpi = 600)
ggsave("plot_N200.pdf", plot = plot_N200, width = 12.8, height = 10.6, units = "in", dpi = 600)

# CP v B ---- ##################################################################

df_CPvB <- sim_data_all %>%
  pivot_longer(
    cols = c(width_exact_CP_mean, width_exact_B_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_CPvB"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "exact_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_CPvB)) %>%
  filter(Se == 0.9 & Sp == 0.9)

method_colors_CPvB <- c("exact_CP" = "#E41A1C",
                        "exact_B" = "#E41A1C")

method_linetypes_CPvB <- c("exact_CP" = "2",
                           "exact_B" = "1")

df_CPvB <- df_CPvB %>%
  mutate(linetype = method_linetypes_CPvB[as.character(method)])

plot_CPvB_50 <- 
  ggplot(df_CPvB %>% filter(N == 50),
         aes(x = true_prev, y = width_CPvB, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  labs(
    x = "True prevalence",
    y = "Interval width",
    title = "Interval width difference between\nClopper-Pearson and Blaker methods",
    subtitle = expression("Population size: 50 (sample size in panels), " * italic("Se") * " = " * italic("Sp") * " = 0.9")
  ) +
  scale_color_manual(values = method_colors_CPvB,
                     labels = c(
                       "Exact Clopper-Pearson CI",
                       "Exact Blaker CI"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(5, 1)), nrow = 1),
         linetype = "none")

ggsave("plot_CPvB.svg", plot = plot_CPvB_50, width = 7, height = 4, units = "in", dpi = 600)
ggsave("plot_CPvB.pdf", plot = plot_CPvB_50, width = 7, height = 4, units = "in", dpi = 600)

# width and point estimate sd ---- #############################################

df_wid_sd <- sim_data_all %>%
  pivot_longer(
    cols = c(width_exact_B_sd, width_simul_B_sd, width_Bayes_sd,
             width_FPC_sd, width_Ge_sd, width_binom_B_sd),
    names_to = "method",
    names_pattern = "width_(.*)_sd",
    values_to = "width_sd"
  ) %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "simul_B",
                                            "Bayes",
                                            "FPC",
                                            "Ge",
                                            "binom_B"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_sd)) %>%
  filter(Se == 0.9 & Sp == 0.9)

df_wid_sd <- df_wid_sd %>%
  mutate(linetype = method_linetypes[as.character(method)])

plot_wid_sd_50 <- 
  ggplot(df_wid_sd %>% filter(N == 50),
         aes(x = true_prev, y = width_sd, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  labs(
    x = "True prevalence",
    y = "SD of CI width"
  ) +
  scale_color_manual(values = method_colors) +
  theme_bw() +
  theme(legend.position = "none")

df_loc_sd <- sim_data_all %>%
  pivot_longer(
    cols = c(loc_exact_sd, loc_simul_sd, loc_Bayes_sd,
             loc_Ge_sd, loc_RG_sd),
    names_to = "method",
    names_pattern = "loc_(.*)_sd",
    values_to = "loc_sd"
  ) %>%
  mutate(method = factor(method, levels = c("exact",
                                            "simul",
                                            "Bayes",
                                            "Ge",
                                            "RG"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, loc_sd)) %>%
  filter(Se == 0.9 & Sp == 0.9)

df_loc_sd <- df_loc_sd %>%
  mutate(linetype = method_linetypes_loc[as.character(method)])

plot_loc_sd_50 <- 
  ggplot(df_loc_sd %>% filter(N == 50),
         aes(x = true_prev, y = loc_sd, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = linetype), position = position_dodge(width = 0.01)) +
  facet_grid(. ~ n, scales = "free_x") +
  labs(
    x = "True prevalence",
    y = "SD of point estimate"
  ) +
  scale_color_manual(values = method_colors_loc) +
  theme_bw() +
  theme(legend.position = "none")

title_grob <- textGrob("Standard deviation of interval width and point estimate\nPopulation size: 50 (sample size in panels)\n",
                       gp = gpar(fontsize = 15, fontface = "bold"))

subheadings <- arrangeGrob(
  textGrob(expression(italic(Se)*" = "*italic(Sp)*" = "*0.9), gp = gpar(fontsize = 12)),
  ncol = 1
)

df_cov_N200_2 <- df_cov_N200 %>%
  mutate(method = factor(method, levels = c("exact_B",
                                            "Bayes",
                                            "simul_B",
                                            "FPC",
                                            "Ge",
                                            "binom_B")))

legend_plot <-
  ggplot(df_cov_N200_2,
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Bayes hypgeom CI",
                       "Simul hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)",
                       "Binom CI (Reiczigel et al. 2010) / Rogan-Gladen estimate"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 2, 2, 2)), nrow = 3))

legend <- cowplot::get_legend(legend_plot)

plot_sd <- 
  grid.arrange(subheadings,
               arrangeGrob(plot_wid_sd_50, plot_loc_sd_50,
                           nrow = 2),
               legend,
               ncol = 1,
               heights = c(0.5, 10, 3),
               top = title_grob)

ggsave("plot_sd.svg", plot = plot_sd, width = 7, height = 7.5, units = "in", dpi = 600)
ggsave("plot_sd.pdf", plot = plot_sd, width = 7, height = 7.5, units = "in", dpi = 600)

# min coverage ---- ############################################################

min(sim_data_all$coverage_exact_CP_mean)
min(sim_data_all$coverage_exact_B_mean)
min(sim_data_all$coverage_simul_CP_mean)
min(sim_data_all$coverage_simul_B_mean)
min(sim_data_all$coverage_Bayes_mean)
min(sim_data_all$coverage_FPC_mean)
min(sim_data_all$coverage_Ge_mean)
min(sim_data_all$coverage_binom_CP_mean)
min(sim_data_all$coverage_binom_B_mean)
