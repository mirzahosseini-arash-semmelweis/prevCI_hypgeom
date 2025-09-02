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

method_colors <- c("exact" = "#E41A1C",
                   "exact_CP" = "#E41A1C",
                   "simul" = "#7CAE00",
                   "simul_CP" = "#7CAE00",
                   "Bayes" = "#00BFC4",
                   "FPC" = "#C77CFF")

sim_data_N50 <- sim_data_all %>% filter(N == 50)

df_cov_N50 <- sim_data_N50 %>%
  pivot_longer(
    cols = c(coverage_exact_CP_mean, coverage_simul_CP_mean, coverage_Bayes_mean,
             coverage_FPC_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

plot_cov_50_09 <-
  ggplot(df_cov_N50 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Coverage",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_50_099 <-
  ggplot(df_cov_N50 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N50 <- sim_data_N50 %>%
  pivot_longer(
    cols = c(width_exact_CP_mean, width_simul_CP_mean, width_Bayes_mean,
             width_FPC_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

plot_wid_50_09 <-
  ggplot(df_wid_N50 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
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
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
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
                       "FPC Wald binom CI"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 2, 1, 2)), nrow = 1))

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
                         ncol = 2),
             legend,
             ncol = 1,
             heights = c(0.5, 10, 1),
             top = title_grob)

ggsave("plot_N50.svg", plot = plot_N50, width = 12.8, height = 6.5, units = "in", dpi = 600)
ggsave("plot_N50.pdf", plot = plot_N50, width = 12.8, height = 6.5, units = "in", dpi = 600)

# N100 ---- ####################################################################

sim_data_N100 <- sim_data_all %>% filter(N == 100)

df_cov_N100 <- sim_data_N100 %>%
  pivot_longer(
    cols = c(coverage_exact_CP_mean, coverage_simul_CP_mean, coverage_Bayes_mean,
             coverage_FPC_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

plot_cov_100_09 <-
  ggplot(df_cov_N100 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Coverage",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_100_099 <-
  ggplot(df_cov_N100 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N100 <- sim_data_N100 %>%
  pivot_longer(
    cols = c(width_exact_CP_mean, width_simul_CP_mean, width_Bayes_mean,
             width_FPC_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

plot_wid_100_09 <-
  ggplot(df_wid_N100 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
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
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
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
                       "FPC Wald binom CI"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 2, 1, 2)), nrow = 1))

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
                           ncol = 2),
               legend,
               ncol = 1,
               heights = c(0.5, 10, 1),
               top = title_grob)

ggsave("plot_N100.svg", plot = plot_N100, width = 12.8, height = 6.5, units = "in", dpi = 600)
ggsave("plot_N100.pdf", plot = plot_N100, width = 12.8, height = 6.5, units = "in", dpi = 600)

# N200 ---- ####################################################################

sim_data_N200 <- sim_data_all %>% filter(N == 200)

df_cov_N200 <- sim_data_N200 %>%
  pivot_longer(
    cols = c(coverage_exact_CP_mean, coverage_simul_CP_mean, coverage_Bayes_mean,
             coverage_FPC_mean),
    names_to = "method",
    names_pattern = "coverage_(.*)_mean",
    values_to = "coverage_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, coverage_mean))

plot_cov_200_09 <-
  ggplot(df_cov_N200 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "Coverage",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

plot_cov_200_099 <-
  ggplot(df_cov_N200 %>%
           filter(Se == 0.99, Sp == 0.99),
         aes(x = true_prev, y = coverage_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0.9, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
  geom_hline(yintercept = 0.95, color = "#B0C4DE", linetype = "longdash", alpha = 0.95) +
  theme_bw() +
  theme(legend.position = "none")

df_wid_N200 <- sim_data_N200 %>%
  pivot_longer(
    cols = c(width_exact_CP_mean, width_simul_CP_mean, width_Bayes_mean,
             width_FPC_mean),
    names_to = "method",
    names_pattern = "width_(.*)_mean",
    values_to = "width_mean"
  ) %>%
  mutate(method = factor(method, levels = c("exact_CP",
                                            "simul_CP",
                                            "Bayes",
                                            "FPC"))) %>%
  select(c(N, true_prev, n, Se, Sp, repl, method, width_mean))

plot_wid_200_09 <-
  ggplot(df_wid_N200 %>%
           filter(Se == 0.9, Sp == 0.9),
         aes(x = true_prev, y = width_mean, color = method, group = interaction(method, n))) +
  geom_point(shape = 1, position = position_dodge(width = 0.01)) +
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
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
  geom_line(aes(linetype = method), position = position_dodge(width = 0.01)) +
  scale_linetype_manual(values = c("exact_CP" = "solid",
                                   "simul_CP" = "dashed",
                                   "Bayes" = "solid",
                                   "FPC" = "dashed")) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(. ~ n, scales = "free_x") +
  scale_color_manual(values = method_colors) +
  labs(y = "",
       x = "True prevalence") +
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
                       "FPC Wald binom CI"
                     ),
                     name = "method") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 2, 1, 2)), nrow = 1))

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
                           ncol = 2),
               legend,
               ncol = 1,
               heights = c(0.5, 10, 1),
               top = title_grob)

ggsave("plot_N200.svg", plot = plot_N200, width = 12.8, height = 6.5, units = "in", dpi = 600)
ggsave("plot_N200.pdf", plot = plot_N200, width = 12.8, height = 6.5, units = "in", dpi = 600)

# CP v B ---- ##################################################################

sim_data_N50[, B_short := (width_exact_CP_mean - width_exact_B_mean)/width_exact_CP_mean]
sim_data_N100[, B_short := (width_exact_CP_mean - width_exact_B_mean)/width_exact_CP_mean]
sim_data_N200[, B_short := (width_exact_CP_mean - width_exact_B_mean)/width_exact_CP_mean]
cell_50_09 <- paste0(
  round(100*min(filter(sim_data_N50, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "-",
  round(100*max(filter(sim_data_N50, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "%"
)
cell_50_099 <- paste0(
  round(100*min(filter(sim_data_N50, Se == 0.99 & Sp == 0.99)$B_short), 1),
  "-",
  round(100*max(filter(sim_data_N50, Se == 0.99 & Sp == 0.99)$B_short), 1),
  "%"
)
cell_100_09 <- paste0(
  round(100*min(filter(sim_data_N100, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "-",
  round(100*max(filter(sim_data_N100, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "%"
)
cell_100_099 <- paste0(
  round(100*min(filter(sim_data_N100, Se == 0.99 & Sp == 0.99)$B_short), 1),
  "-",
  round(100*max(filter(sim_data_N100, Se == 0.99 & Sp == 0.99)$B_short), 1),
  "%"
)
cell_200_09 <- paste0(
  round(100*min(filter(sim_data_N200, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "-",
  round(100*max(filter(sim_data_N200, Se == 0.9 & Sp == 0.9)$B_short), 1),
  "%"
)
cell_200_099 <- paste0(
  format(round(100*min(filter(sim_data_N200, Se == 0.99 & Sp == 0.99)$B_short), 1), nsmall = 1),
  "-",
  round(100*max(filter(sim_data_N200, Se == 0.99 & Sp == 0.99)$B_short), 1),
  "%"
)
CPvB_tab <- data.frame(
                       "50" = c(cell_50_09, cell_50_099),
                       "100" = c(cell_100_09, cell_100_099),
                       "200" = c(cell_200_09, cell_200_099),
                       row.names = c("Se = Sp = 0.9", "Se = Sp = 0.99")
                       )
colnames(CPvB_tab) <- c("N = 50", "N = 100", "N = 200")
write.csv(CPvB_tab, "CPvB_tab.csv", row.names = TRUE)

# point estimate bias and sd ---- ##############################################

sim_data_N50_S09 <- subset(sim_data_N50, Se == 0.9 & Sp == 0.9 & true_prev %in% c(0.01, 0.1, 0.5))
col_N50_S09 <- paste0("$",
                      round(sim_data_N50_S09$loc_exact_mean - sim_data_N50_S09$true_prev, 3),
                      " / ",
                      round(sim_data_N50_S09$loc_exact_sd, 3),
                      "$")
sim_data_N50_S099 <- subset(sim_data_N50, Se == 0.99 & Sp == 0.99 & true_prev %in% c(0.01, 0.1, 0.5))
col_N50_S099 <- paste0("$",
                       round(sim_data_N50_S099$loc_exact_mean - sim_data_N50_S099$true_prev, 3),
                       " / ",
                       round(sim_data_N50_S099$loc_exact_sd, 3),
                       "$")
sim_data_N100_S09 <- subset(sim_data_N100, Se == 0.9 & Sp == 0.9 & true_prev %in% c(0.01, 0.1, 0.5))
col_N100_S09 <- paste0("$",
                       round(sim_data_N100_S09$loc_exact_mean - sim_data_N100_S09$true_prev, 3),
                       " / ",
                       round(sim_data_N100_S09$loc_exact_sd, 3),
                       "$")
sim_data_N100_S099 <- subset(sim_data_N100, Se == 0.99 & Sp == 0.99 & true_prev %in% c(0.01, 0.1, 0.5))
col_N100_S099 <- paste0("$",
                        round(sim_data_N100_S099$loc_exact_mean - sim_data_N100_S099$true_prev, 3),
                        " / ",
                        round(sim_data_N100_S099$loc_exact_sd, 3),
                        "$")
sim_data_N200_S09 <- subset(sim_data_N200, Se == 0.9 & Sp == 0.9 & true_prev %in% c(0.01, 0.1, 0.5))
col_N200_S09 <- paste0("$",
                       round(sim_data_N200_S09$loc_exact_mean - sim_data_N200_S09$true_prev, 3),
                       " / ",
                       round(sim_data_N200_S09$loc_exact_sd, 3),
                       "$")
sim_data_N200_S099 <- subset(sim_data_N200, Se == 0.99 & Sp == 0.99 & true_prev %in% c(0.01, 0.1, 0.5))
col_N200_S099 <- paste0("$",
                        round(sim_data_N200_S099$loc_exact_mean - sim_data_N200_S099$true_prev, 3),
                        " / ",
                        round(sim_data_N200_S099$loc_exact_sd, 3),
                        "$")
loc_tab <- data.frame(
                      "prev" = rep(c(0.01, 0.1, 0.5), 3),
                      "n%" = rep(c("10%", "20%", "50%"), each = 3),
                      "col_N50_S09" = col_N50_S09,
                      "col_N50_S099" = col_N50_S099,
                      "col_N100_S09" = col_N100_S09,
                      "col_N100_S099" = col_N100_S099,
                      "col_N200_S09" = col_N200_S09,
                      "col_N200_S099" = col_N200_S099
                      )
write.csv(loc_tab, "loc_tab.csv", row.names = FALSE, quote = TRUE)