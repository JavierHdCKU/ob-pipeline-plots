#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(stringr)
  library(viridis)
  library(scales)
  library(jsonlite)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--f1_cv"), type = "character", help = "f1_macro_by_crossvalidation.tsv"),
  make_option(c("--f1_weighted"), type = "character", help = "samples_vs_f1_weighted.tsv"),
  make_option(c("--conf_matrix"), type = "character", help = "per_population_confusion.tsv"),
  make_option(c("--perf_tsv"), type = "character", help = "performances.tsv"),
  make_option(c("--meta_json"), type = "character", help = "dataset_metadata.json"),
  make_option(c("-o", "--out_dir"), type = "character", default = "./plots", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Ensure output directory exists
if (!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# SHARED MAPPINGS & THEMES
# ------------------------------------------------------------------------------
name_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 'Chikungunya virus infection',
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 'Healthy',
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 'COVID-19 infection',
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 'Stim Blood',
  'dataset_name-Samusik_seed-42'                    = 'Mouse Bone Marrow',
  'dataset_name-Transformed_seed-42'                = 'Transformed',
  'dataset_name-flowcyt_seed-42' = 'Human Bone Marrow',
  'dataset_name-Levine_seed-42' = 'Levine'
)

marker_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = 37,
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = 24,
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = 24,
  'dataset_name-FR-FCM-Z3YR_seed-42'                = 38,
  'dataset_name-Samusik_seed-42'                    = 39,
  'dataset_name-Transformed_seed-42'                = 33,
  'dataset_name-flowcyt_seed-42' = 12,
  'dataset_name-Levine_seed-42' = 32
)

platform_map <- c(
  'dataset_name-FR-FCM-Z238_infection_final_seed-42' = "CyTOF",
  'dataset_name-FR-FCM-Z2KP_healthy_final_seed-42'  = "FlowCyt",
  'dataset_name-FR-FCM-Z2KP_virus_final_seed-42'    = "FlowCyt",
  'dataset_name-FR-FCM-Z3YR_seed-42'                = "CyTOF",
  'dataset_name-Samusik_seed-42'                    = "CyTOF",
  'dataset_name-Transformed_seed-42'                = "CyTOF",
  'dataset_name-flowcyt_seed-42' = "FlowCyt",
  'dataset_name-Levine_seed-42' = 'CyTOF'
)

# ------------------------------------------------------------------------------
# FIGURE 1: HEATMAP & COMPOSITE BOXPLOTS
# ------------------------------------------------------------------------------
cat("Generating Figure 1 (Composite Heatmap)...\n")
df1 <- read_tsv(opt$f1_cv, show_col_types = FALSE)

df1_clean <- df1 %>%
  filter(!str_detect(dataset, regex("sub-sampling", ignore_case = TRUE))) %>%
  mutate(
    platform = platform_map[dataset],
    markers  = marker_map[dataset],
    clean_name = recode(dataset, !!!name_map),
    display_name = ifelse(!is.na(markers), 
                          paste0(clean_name, " (", markers, ")"), 
                          clean_name),
    platform = ifelse(is.na(platform), "Other", platform)
  )

df1_no_random <- df1_clean %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

# 3. DYNAMIC ORDERING
dataset_order_df1 <- df1_no_random %>%
  group_by(platform, display_name) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE), .groups = "drop") %>%
  arrange(platform, global_mean)

dataset1_order <- dataset_order_df1$display_name

tool_order1 <- df1_clean %>%
  group_by(model) %>%
  summarise(global_mean = mean(f1_macro, na.rm = TRUE), .groups = "drop") %>%
  arrange(global_mean) %>%
  pull(model)

df1_clean$model <- factor(df1_clean$model, levels = tool_order1)
df1_clean$display_name <- factor(df1_clean$display_name, levels = dataset1_order)
df1_no_random$display_name <- factor(df1_no_random$display_name, levels = dataset1_order)

heatmap_data1 <- df1_clean %>%
  group_by(model, display_name, platform) %>%
  summarise(mean_f1 = mean(f1_macro, na.rm = TRUE), .groups = "drop")

# 4. PLOTTING COMPONENTS
theme_gb_base1 <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    plot.margin = margin(2, 2, 2, 2)
  )

p1_heatmap <- ggplot(heatmap_data1, aes(x = display_name, y = model, fill = mean_f1)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_f1)), size = 1.8, color = "black") +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "viridis", limits = c(0, 1), guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Training set", y = NULL) +
  theme_gb_base1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.title.x = element_text(margin = margin(b = 20), size = 9), 
    panel.spacing = unit(8, "pt"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid = element_blank()
  )

p1_right <- ggplot(df1_clean, aes(x = f1_macro, y = model)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "F1-score", y = NULL) +
  theme_gb_base1 +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.y = element_blank(),
    plot.margin = margin(l = 2)
  )

p1_bottom <- ggplot(df1_no_random, aes(x = display_name, y = f1_macro)) +
  geom_boxplot(outlier.size = 0.1, lwd = 0.25, fill = "white") +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "F1-score") +
  theme_gb_base1 +
  theme(
    panel.background = element_rect(fill = "grey95", color = NA),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(6, "pt"),
    plot.margin = margin(t = 2)
  )

# 5. COMPOSITION & SAVING
design1 <- "
AAAAB
AAAAB
AAAAB
CCCC#
"
final_plot1 <- wrap_plots(A = p1_heatmap, B = p1_right, C = p1_bottom, design = design1)

n_tools1 <- length(unique(df1_clean$model))
n_datasets1 <- length(unique(df1_clean$display_name))
final_h1 <- min(max((n_tools1 * 8) + 80, 120), 225)
final_w1 <- min(max((n_datasets1 * 15) + 80, 140), 180)

output_fig1 <- file.path(opt$out_dir, "Figure1_Heatmap_Composite.png")
ggsave(output_fig1, plot = final_plot1, width = final_w1, height = final_h1, units = "mm", dpi = 600)

cat(paste0("Success: Figure 1 saved to ", output_fig1, "\n"))

# ------------------------------------------------------------------------------
# FIGURE 2: PER-POPULATION BOXPLOTS (DISTRIBUTION ACROSS CV RUNS)
# ------------------------------------------------------------------------------
cat("Generating Figure 2 (Individual Dataset Plots)...\n")
df2 <- read_tsv(opt$conf_matrix, show_col_types = FALSE)

df2_clean <- df2 %>%
  mutate(dataset = recode(dataset, !!!name_map)) %>%
  mutate(
    tp = as.numeric(tp),
    fp = as.numeric(fp),
    fn = as.numeric(fn),
    denom = (2 * tp) + fp + fn,
    f1_score = ifelse(denom > 0, (2 * tp) / denom, 0)
  )

# Define Style for Figure 2
theme_gb_white2 <- theme_classic(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(4, "mm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

unique_datasets2 <- unique(df2_clean$dataset)
cat(paste("Found", length(unique_datasets2), "datasets. Generating plots...\n"))

for (ds_name in unique_datasets2) {
  
  plot_data2 <- df2_clean %>% filter(dataset == ds_name)
  n_pops2 <- length(unique(plot_data2$population_name))
  n_tools2 <- length(unique(plot_data2$model))
  
  p2 <- ggplot(plot_data2, aes(x = population_name, y = f1_score, fill = model)) +
    geom_boxplot(
      position = position_dodge(width = 0.8), 
      outlier.size = 0.3, 
      lwd = 0.2, 
      color = "black",
      alpha = 0.9
    ) +
    scale_fill_viridis_d(option = "turbo", name = "Classifiers") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = ds_name, x = "Cell Populations", y = "F1 Score") +
    theme_gb_white2 +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Dynamic Width Calculation based on complexity
  required_width_mm2 <- (n_pops2 * n_tools2 * 5) + 50
  final_width2 <- min(max(required_width_mm2, 100), 400)
  
  safe_name2 <- str_replace_all(ds_name, "[^A-Za-z0-9]", "_")
  filename2 <- file.path(opt$out_dir, paste0("Figure2_", safe_name2, ".png"))
  
  ggsave(filename2, plot = p2, width = final_width2, height = 100, 
         units = "mm", dpi = 600, device = "png")
  
  cat(paste("Saved:", filename2, "\n"))
}

cat("Success: All Figure 2 plots have been generated.\n")

# ------------------------------------------------------------------------------
# SHARED METADATA EXTRACTION & THEME FOR FIG 3 & 4
# ------------------------------------------------------------------------------
cat("Processing Metadata for Figures 3 & 4...\n")
metadata_json <- fromJSON(opt$meta_json)

extract_metadata_stats <- function(meta_list) {
  ds_ids <- names(meta_list)
  ds_ids_filtered <- ds_ids[!str_detect(ds_ids, regex("Levine", ignore_case = TRUE))]
  
  results <- lapply(ds_ids_filtered, function(id) {
    entry <- meta_list[[id]]
    cells <- as.numeric(entry$cells_per_sample)
    avg_cells <- if(all(is.na(cells))) NA else mean(cells, na.rm = TRUE)
    
    data.frame(
      dataset_id    = id,
      n_markers     = if(id %in% names(marker_map)) as.numeric(marker_map[id]) else NA,
      n_samples     = as.numeric(entry$sample_count),
      n_populations = as.numeric(entry$population_count),
      mean_cells    = avg_cells,
      stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, results))
}

df_metadata_shared <- extract_metadata_stats(metadata_json)

theme_gb_scatter_shared <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# FIGURE 3: PERFORMANCE vs METADATA (ACCURACY)
# ------------------------------------------------------------------------------
cat("Generating Figure 3...\n")
df3_raw <- read_tsv(opt$f1_weighted, show_col_types = FALSE)

df3_avg <- df3_raw %>%
  group_by(dataset, model) %>%
  summarise(mean_f1_weighted = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop")

df3_plot <- df3_avg %>%
  filter(!str_detect(dataset, regex("Levine", ignore_case = TRUE))) %>%
  left_join(df_metadata_shared, by = c("dataset" = "dataset_id")) %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

create_fig3_panel <- function(data, x_var, x_label, is_log = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_f1_weighted, color = model)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.6) +
    geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2, alpha = 0.8) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = x_label, y = "Mean F1 (Weighted)") +
    theme_gb_scatter_shared
  
  if(is_log) {
    p <- p + scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x))
  }
  return(p)
}

p3_cells   <- create_fig3_panel(df3_plot, "mean_cells", "Mean Cells / Sample", is_log = TRUE)
p3_markers <- create_fig3_panel(df3_plot, "n_markers", "Number of Markers")
p3_samples <- create_fig3_panel(df3_plot, "n_samples", "Number of Samples")
p3_pops    <- create_fig3_panel(df3_plot, "n_populations", "Number of Populations")

final_fig3 <- (p3_cells + p3_markers) / (p3_samples + p3_pops) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"), plot.tag = element_text(face = "bold"))

ggsave(file.path(opt$out_dir, "Figure3_Accuracy_Metadata.png"), final_fig3, width = 180, height = 160, units = "mm", dpi = 600)

# ------------------------------------------------------------------------------
# FIGURE 4: PERFORMANCE vs METADATA (RUNTIME)
# ------------------------------------------------------------------------------
cat("Generating Figure 4...\n")
df4_perf <- read_tsv(opt$perf_tsv, show_col_types = FALSE)

df4_time_avg <- df4_perf %>%
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  mutate(
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[,2],
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[,2],
    dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed)
  ) %>%
  filter(!is.na(extracted_name)) %>%
  group_by(dataset_id, model = module) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

df4_plot <- df4_time_avg %>%
  inner_join(df_metadata_shared, by = "dataset_id") %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

create_fig4_panel <- function(data, x_var, x_label, is_log_x = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_time_sec, color = model)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.6) +
    geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2, alpha = 0.8) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
    labs(x = x_label, y = "Mean Runtime (s)") +
    theme_gb_scatter_shared
  
  if(is_log_x) {
    p <- p + scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x))
  }
  return(p)
}

p4_cells   <- create_fig4_panel(df4_plot, "mean_cells", "Mean Cells / Sample", is_log_x = TRUE)
p4_markers <- create_fig4_panel(df4_plot, "n_markers", "Number of Markers")
p4_samples <- create_fig4_panel(df4_plot, "n_samples", "Number of Samples")
p4_pops    <- create_fig4_panel(df4_plot, "n_populations", "Number of Populations")

final_fig4 <- (p4_cells + p4_markers) / (p4_samples + p4_pops) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"), plot.tag = element_text(face = "bold"))

ggsave(file.path(opt$out_dir, "Figure4_Runtime_Metadata.png"), final_fig4, width = 180, height = 160, units = "mm", dpi = 600)

# ------------------------------------------------------------------------------
# FIGURE 5: SCALABILITY (SUB-SAMPLING)
# ------------------------------------------------------------------------------
cat("Generating Figure 5 (Scalability Panels)...\n")

# 1. Helper for Dataset String Parsing
extract_size <- function(name) {
  as.numeric(str_extract(name, "(?<=sub-sampling-)\\d+"))
}

# 2. Process Performance (Time)
df5_time <- read_tsv(opt$perf_tsv, show_col_types = FALSE) %>%
  filter(str_detect(params, "sub-sampling")) %>%
  mutate(
    train_size = as.numeric(str_match(params, 'sub-sampling[^:]+:\\s*\"+([^\"]+)')[,2]),
    model = module
  ) %>%
  filter(!model %in% c("metrics", "flow_metrics", "data_import", "data_preprocessing", "random")) %>%
  group_by(model, train_size) %>%
  summarise(mean_time = mean(s, na.rm = TRUE), .groups = "drop")

# 3. Process F1 (Accuracy)
df5_f1 <- read_tsv(opt$f1_weighted, show_col_types = FALSE) %>%
  filter(str_detect(dataset, "sub-sampling")) %>%
  mutate(train_size = extract_size(dataset)) %>%
  group_by(model, train_size) %>%
  summarise(mean_f1 = mean(f1_weighted_median, na.rm = TRUE), .groups = "drop") %>%
  filter(!str_detect(model, "random"))

# 4. Harmonize Colors across both dataframes
df5_validation <- full_join(df5_time, df5_f1, by = c("model", "train_size")) %>%
  arrange(model, train_size)

my_colors5 <- setNames(viridis::turbo(length(unique(df5_validation$model))), 
                       unique(df5_validation$model))

theme_gb_scalability <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5)
  )

# 5. Plotting
p5_time <- ggplot(df5_time, aes(x = train_size, y = mean_time, color = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(shape = 21, fill = "white", stroke = 1, size = 2) +
  scale_y_log10(labels = label_log()) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  scale_color_manual(values = my_colors5, name = "Classifiers") +
  labs(title = "Computational Scalability", x = "Training Set Size", y = "Runtime (s)") +
  theme_gb_scalability +
  theme(legend.position = "none")

p5_f1 <- ggplot(df5_f1, aes(x = train_size, y = mean_f1, color = model, group = model)) +
  geom_line(linewidth = 0.7) +
  geom_point(shape = 21, fill = "white", stroke = 1, size = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(labels = label_number(suffix = "K", scale = 1e-3), breaks = c(10000, 50000, 100000)) +
  scale_color_manual(values = my_colors5, name = "Classifiers") +
  labs(title = "Classification Stability", x = "Training Set Size", y = "Mean F1 (Weighted)") +
  theme_gb_scatter_shared # Reusing the shared scatter theme for consistency

# 6. Combine & Save
final_fig5 <- (p5_time / p5_f1) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "bottom")

output_fig5 <- file.path(opt$out_dir, "Figure5_Scalability.png")
ggsave(output_fig5, final_fig5, width = 120, height = 180, units = "mm", dpi = 600)

cat(paste0("Success: Figure 5 saved to ", output_fig5, "\n"))
cat(paste0("\n[ALL TASKS COMPLETE] All plots saved to: ", opt$out_dir, "\n"))
###
# Usage example:
#Rscript generate_all_figures.r \
#  --f1_cv /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/f1_macro_by_crossvalidation.tsv \
#  --f1_weighted /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/samples_vs_f1_weighted.tsv \
#  --conf_matrix /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv \
#  --perf_tsv /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/performances.tsv \
#  --meta_json /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --out_dir /home/projects/dp_immunoth/people/javher/projects/benchmarking/ob-pipeline-plots
###