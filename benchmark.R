# benchmarking
library(microbenchmark)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)

bench <- microbenchmark(
  exact50 = prevCI.hypgeom.exact(50, 5, 3, 0.9, 0.9),
  exact75 = prevCI.hypgeom.exact(75, 8, 4, 0.9, 0.9),
  exact100 = prevCI.hypgeom.exact(100, 10, 5, 0.9, 0.9),
  exact125 = prevCI.hypgeom.exact(125, 13, 7, 0.9, 0.9),
  exact150 = prevCI.hypgeom.exact(150, 15, 8, 0.9, 0.9),
  exact200 = prevCI.hypgeom.exact(200, 20, 10, 0.9, 0.9),
  exact250 = prevCI.hypgeom.exact(250, 25, 13, 0.9, 0.9),
  exact300 = prevCI.hypgeom.exact(300, 30, 15, 0.9, 0.9),
  simul50 = prevCI.hypgeom.simul(50, 5, 3, 0.9, 0.9),
  simul75 = prevCI.hypgeom.simul(75, 8, 4, 0.9, 0.9),
  simul100 = prevCI.hypgeom.simul(100, 10, 5, 0.9, 0.9),
  simul125 = prevCI.hypgeom.simul(125, 13, 7, 0.9, 0.9),
  simul150 = prevCI.hypgeom.simul(150, 15, 8, 0.9, 0.9),
  simul200 = prevCI.hypgeom.simul(200, 20, 10, 0.9, 0.9),
  simul250 = prevCI.hypgeom.simul(250, 25, 13, 0.9, 0.9),
  simul300 = prevCI.hypgeom.simul(300, 30, 15, 0.9, 0.9),
  Bayes50 = prevCI.hypgeom.Bayes(50, 5, 3, 0.9, 0.9),
  Bayes75 = prevCI.hypgeom.Bayes(75, 8, 4, 0.9, 0.9),
  Bayes100 = prevCI.hypgeom.Bayes(100, 10, 5, 0.9, 0.9),
  Bayes125 = prevCI.hypgeom.Bayes(125, 13, 7, 0.9, 0.9),
  Bayes150 = prevCI.hypgeom.Bayes(150, 15, 8, 0.9, 0.9),
  Bayes200 = prevCI.hypgeom.Bayes(200, 20, 10, 0.9, 0.9),
  Bayes250 = prevCI.hypgeom.Bayes(250, 25, 13, 0.9, 0.9),
  Bayes300 = prevCI.hypgeom.Bayes(300, 30, 15, 0.9, 0.9),
  Ge50 = RS_BC(50, 5, 3, 0.9, 0.9),
  Ge75 = RS_BC(75, 8, 4, 0.9, 0.9),
  Ge100 = RS_BC(100, 10, 5, 0.9, 0.9),
  Ge125 = RS_BC(125, 13, 7, 0.9, 0.9),
  Ge150 = RS_BC(150, 15, 8, 0.9, 0.9),
  Ge200 = RS_BC(200, 20, 10, 0.9, 0.9),
  Ge250 = RS_BC(250, 25, 13, 0.9, 0.9),
  Ge300 = RS_BC(300, 30, 15, 0.9, 0.9),
  FPC50 = RS_Wald(50, 5, 3, 0.9, 0.9),
  FPC75 = RS_Wald(75, 8, 4, 0.9, 0.9),
  FPC100 = RS_Wald(100, 10, 5, 0.9, 0.9),
  FPC125 = RS_Wald(125, 13, 7, 0.9, 0.9),
  FPC150 = RS_Wald(150, 15, 8, 0.9, 0.9),
  FPC200 = RS_Wald(200, 20, 10, 0.9, 0.9),
  FPC250 = RS_Wald(250, 25, 13, 0.9, 0.9),
  FPC300 = RS_Wald(300, 30, 15, 0.9, 0.9),
  times = 1
)

benchmark <- data.frame(method = str_extract(bench$expr, "^[A-Za-z]+"),
                        N = as.integer(str_extract(bench$expr, "[0-9]+")),
                        time = bench$time)
benchmark <- benchmark %>%
  mutate(method = recode(method,
                         "exact" = "exact",
                         "simul" = "simul",
                         "Bayes" = "Bayes",
                         "FPC" = "FPC",
                         "Ge" = "Ge"),
         method = factor(method, levels = c("exact", "Bayes", "simul", "FPC", "Ge")))

bench2 <- microbenchmark(
  exact50 = prevCI.hypgeom.exact(50, 25, 13, 0.9, 0.9),
  exact75 = prevCI.hypgeom.exact(75, 38, 19, 0.9, 0.9),
  exact100 = prevCI.hypgeom.exact(100, 50, 25, 0.9, 0.9),
  exact125 = prevCI.hypgeom.exact(125, 63, 32, 0.9, 0.9),
  exact150 = prevCI.hypgeom.exact(150, 75, 38, 0.9, 0.9),
  exact200 = prevCI.hypgeom.exact(200, 100, 50, 0.9, 0.9),
  exact250 = prevCI.hypgeom.exact(250, 125, 63, 0.9, 0.9),
  exact300 = prevCI.hypgeom.exact(300, 150, 75, 0.9, 0.9),
  simul50 = prevCI.hypgeom.simul(50, 25, 13, 0.9, 0.9),
  simul75 = prevCI.hypgeom.simul(75, 38, 19, 0.9, 0.9),
  simul100 = prevCI.hypgeom.simul(100, 50, 25, 0.9, 0.9),
  simul125 = prevCI.hypgeom.simul(125, 63, 32, 0.9, 0.9),
  simul150 = prevCI.hypgeom.simul(150, 75, 38, 0.9, 0.9),
  simul200 = prevCI.hypgeom.simul(200, 100, 50, 0.9, 0.9),
  simul250 = prevCI.hypgeom.simul(250, 125, 63, 0.9, 0.9),
  simul300 = prevCI.hypgeom.simul(300, 150, 75, 0.9, 0.9),
  Bayes50 = prevCI.hypgeom.Bayes(50, 25, 13, 0.9, 0.9),
  Bayes75 = prevCI.hypgeom.Bayes(75, 38, 19, 0.9, 0.9),
  Bayes100 = prevCI.hypgeom.Bayes(100, 50, 25, 0.9, 0.9),
  Bayes125 = prevCI.hypgeom.Bayes(125, 63, 32, 0.9, 0.9),
  Bayes150 = prevCI.hypgeom.Bayes(150, 75, 38, 0.9, 0.9),
  Bayes200 = prevCI.hypgeom.Bayes(200, 100, 50, 0.9, 0.9),
  Bayes250 = prevCI.hypgeom.Bayes(250, 125, 63, 0.9, 0.9),
  Bayes300 = prevCI.hypgeom.Bayes(300, 150, 75, 0.9, 0.9),
  Ge50 = RS_BC(50, 25, 13, 0.9, 0.9),
  Ge75 = RS_BC(75, 38, 19, 0.9, 0.9),
  Ge100 = RS_BC(100, 50, 25, 0.9, 0.9),
  Ge125 = RS_BC(125, 63, 32, 0.9, 0.9),
  Ge150 = RS_BC(150, 75, 38, 0.9, 0.9),
  Ge200 = RS_BC(200, 100, 50, 0.9, 0.9),
  Ge250 = RS_BC(250, 125, 63, 0.9, 0.9),
  Ge300 = RS_BC(300, 150, 75, 0.9, 0.9),
  FPC50 = RS_BC(50, 25, 13, 0.9, 0.9),
  FPC75 = RS_BC(75, 38, 19, 0.9, 0.9),
  FPC100 = RS_BC(100, 50, 25, 0.9, 0.9),
  FPC125 = RS_BC(125, 63, 32, 0.9, 0.9),
  FPC150 = RS_BC(150, 75, 38, 0.9, 0.9),
  FPC200 = RS_BC(200, 100, 50, 0.9, 0.9),
  FPC250 = RS_BC(250, 125, 63, 0.9, 0.9),
  FPC300 = RS_BC(300, 150, 75, 0.9, 0.9),
  times = 1
)

benchmark2 <- data.frame(method = str_extract(bench2$expr, "^[A-Za-z]+"),
                         N = as.integer(str_extract(bench2$expr, "[0-9]+")),
                         time = bench2$time)
benchmark2 <- benchmark2 %>%
  mutate(method = recode(method,
                         "exact" = "exact",
                         "simul" = "simul",
                         "Bayes" = "Bayes",
                         "FPC" = "FPC",
                         "Ge" = "Ge"),
         method = factor(method, levels = c("exact", "Bayes", "simul", "FPC", "Ge")))

benchmark$pn <- "10%"
benchmark2$pn <- "50%"
benchmark_all <- rbind(benchmark, benchmark2)
benchmark_all$pn <- as.factor(benchmark_all$pn)

metric_colors <- c("exact" = "#E41A1C",
                   "simul" = "#4CBB17",
                   "Bayes" = "#377EB8",
                   "FPC" = "#E41A1C",
                   "Ge" = "#377EB8")

metric_linetypes <- c("exact" = "1",
                      "simul" = "1",
                      "Bayes" = "1",
                      "FPC" = "2",
                      "Ge" = "2")

benchmark_all <- benchmark_all %>%
  mutate(linetype = metric_linetypes[as.character(method)])

bench_plot <-
  ggplot(benchmark_all, aes(x = N, y = time, color = method)) +
  geom_point(shape = 1) +
  geom_line(aes(group = method, linetype = linetype), show.legend = FALSE) +
  scale_color_manual(values = metric_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Simul hypgeom CI",
                       "Bayes hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)"
                     )) +
  facet_grid(. ~ pn) +
  labs(
    title = "Benchmarking of method runtimes",
    subtitle = expression(paste(italic(Se), " = ", italic(Sp), " = 0.9",
                                ", Sample size kept constant at 10% or 50% of population")),
    x = "Population size",
    y = "Time on log scale",
    color = "method"
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none")

legend_plot <-
  ggplot(benchmark_all, aes(x = N, y = time, color = method)) +
  geom_point(shape = 1) +
  geom_line(aes(group = method)) +
  scale_color_manual(values = metric_colors,
                     labels = c(
                       "Exact hypgeom CI",
                       "Bayes hypgeom CI",
                       "Simul hypgeom CI",
                       "FPC Wald binom CI (Ge et al. 2024)",
                       "Bayes CI (Ge et al. 2024)"
                     )) +
  facet_grid(. ~ pn) +
  labs(
    title = "Benchmarking of method runtimes",
    subtitle = expression(paste(italic(Se), " = ", italic(Sp), " = 0.9",
                                ", Sample size kept constant at 10% or 50% of population")),
    x = "Population size",
    y = "Time on log scale",
    color = "method"
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1, 1, 2, 2)), nrow = 3))


legend <- cowplot::get_legend(legend_plot)

plotbench <- 
  grid.arrange(arrangeGrob(bench_plot),
               legend,
               ncol = 1,
               heights = c(10, 3))

ggsave("plotbench.svg", plot = plotbench, width = 8.2, height = 6, units = "in", dpi = 600)
ggsave("plotbench.pdf", plot = plotbench, width = 8.2, height = 6, units = "in", dpi = 600)