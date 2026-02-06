#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(viridis)
  library(scales)
  library(patchwork)
  library(stringr)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--f1_input"), type = "character", default = NULL, 
              help = "Path to samples_vs_f1_weighted.tsv", metavar = "character"),
  make_option(c("-p", "--perf_input"), type = "character", default = NULL, 
              help = "Path to performances.tsv", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure5_Scalability.png", 
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validation
if (is.null(opt$f1_input) || is.null(opt$perf_input)) {
  print_help(opt_parser)
  stop("Both --f1_input and --perf_input must be provided.", call. = FALSE)
}

path_f1 <- opt$f1_input
path_perf <- opt$perf_input
output_file <- opt$output

# ------------------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------------------
df_f1_raw <- read_tsv(path_f1, show_col_types = FALSE)
df_perf_raw <- read_tsv(path_perf, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. DATA PROCESSING (LEVINE ONLY)
# ------------------------------------------------------------------------------

# Helper to extract size from dataset strings
extract_size <- function(name) {
  as.numeric(str_extract(name, "(?<=sub-sampling-)\\d+"))
}

# A. Process Performance (Time)
df_levine_time <- df_perf_raw %>%
  filter(str_detect(params, "sub-sampling")) %>%
  mutate(
    train_size = as.numeric(str_match(params, 'sub-sampling[^:]+:\\s*\"+([^\"]+)')[,2]),
    model = module
  ) %>%
  filter(!model %in% c("metrics", "flow_metrics", "data_import", "data_preprocessing", "random")) %>%
  group_by(model, train_size) %>%
  summarise(mean_time = mean(s, na.rm = TRUE), .groups = "drop")

# B. Process F1 (Accuracy)
df_levine_f1 <- df_f1_raw %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(train_size = extract_size(dataset)) %>%
  group_by(model, train_size) %>%
  summarise(mean_f1 = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop") %>%
  filter(!str_detect(model, "random"))

# Join for color consistency
df_validation <- full_join(df_levine_time, df_levine_f1, by = c("model", "train_size")) %>%
  arrange(model, train_size)

# ------------------------------------------------------------------------------
# 3. HARMONIZED PLOTTING
# ------------------------------------------------------------------------------
my_colors <- setNames(viridis::turbo(length(unique(df_validation$model))), 
                      unique(df_validation$model))

theme_gb_scalability <- theme_minimal(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

p_time <- ggplot(df_levine_time, aes(x = train_size, y = mean_time, color = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(shape = 21, fill = "white", stroke = 1, size = 2) +
  scale_y_log10(labels = label_log()) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  scale_color_manual(values = my_colors, name = "Classifiers") +
  labs(title = "Computational Scalability (Healthy FlowCyt dataset)", x = "Training Set Size", y = "Runtime (s)") +
  theme_gb_scalability +
  theme(legend.position = "none")

p_f1 <- ggplot(df_levine_f1, aes(x = train_size, y = mean_f1, color = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(shape = 21, fill = "white", stroke = 1, size = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  scale_color_manual(values = my_colors, name = "Classifiers") +
  labs(title = "Classification Stability (Healthy FlowCyt dataset)", x = "Training Set Size", y = "Mean F1 (Weighted)") +
  theme_gb_scalability

# ------------------------------------------------------------------------------
# 4. COMBINE & SAVE
# ------------------------------------------------------------------------------
final_fig5 <- (p_time / p_f1) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "bottom")

ggsave(output_file, final_fig5, width = 120, height = 180, units = "mm", dpi = 600)

cat(paste0("Success: Figure 5 saved to ", output_file, "\n"))

###
# Usage example:
#Rscript Fig5_plot.r \
#  --f1_input ./ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#  --perf_input ./ob-blob-metrics/out/performances.tsv \
#  --output ./ob-pipeline-plots/Figure5_Scalability.png
###