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
  library(jsonlite)
  library(stringr)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--perf_input"), type = "character", default = NULL, 
              help = "Path to the performances.tsv file", metavar = "character"),
  make_option(c("-j", "--meta_json"), type = "character", default = NULL, 
              help = "Path to dataset_metadata.json", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "Figure4_Runtime_Metadata.png", 
              help = "Output filename", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validation
if (is.null(opt$perf_input) || is.null(opt$meta_json)) {
  print_help(opt_parser)
  stop("Both --perf_input and --meta_json must be provided.", call. = FALSE)
}

path_perf <- opt$perf_input
path_meta_json <- opt$meta_json
output_file <- opt$output

# ------------------------------------------------------------------------------
# 1. DATA LOADING
# ------------------------------------------------------------------------------
df_perf <- read_tsv(path_perf, show_col_types = FALSE)
metadata_json <- fromJSON(path_meta_json)

# ------------------------------------------------------------------------------
# 2. PARSE PERFORMANCE & SUMMARIZE
# ------------------------------------------------------------------------------
df_time_avg <- df_perf %>%
  filter(!module %in% c("data_import", "data_preprocessing", "flow_metrics", "f1_score", "metric_collector", "metrics")) %>%
  mutate(
    extracted_name = str_match(params, 'dataset_name[^:]+:\\s*\"+([^\"]+)')[,2],
    extracted_seed = str_match(params, 'seed[^:]+:\\s*\"+([^\"]+)')[,2],
    dataset_id = paste0("dataset_name-", extracted_name, "_seed-", extracted_seed)
  ) %>%
  filter(!is.na(extracted_name)) %>%
  group_by(dataset_id, model = module) %>%
  summarise(mean_time_sec = mean(s, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. METADATA EXTRACTION
# ------------------------------------------------------------------------------
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

df_metadata <- extract_metadata_stats(metadata_json)

df_plot <- df_time_avg %>%
  inner_join(df_metadata, by = "dataset_id") %>%
  filter(!str_detect(model, regex("random", ignore_case = TRUE)))

# ------------------------------------------------------------------------------
# 4. THEME & PLOTTING
# ------------------------------------------------------------------------------
theme_gb_scatter <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold", color = "black"),
    panel.background = element_rect(fill = "grey98", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

create_time_plot <- function(data, x_var, x_label, is_log_x = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = mean_time_sec, color = model)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, alpha = 0.6) +
    geom_point(aes(fill = model), shape = 21, color = "black", stroke = 0.2, size = 2, alpha = 0.8) +
    scale_color_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo") +
    scale_y_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x)) +
    labs(x = x_label, y = "Mean Runtime (s)") +
    theme_gb_scatter
  
  if(is_log_x) {
    p <- p + scale_x_log10(labels = label_log(), breaks = trans_breaks("log10", function(x) 10^x))
  }
  return(p)
}

# Generate panels
p_cells   <- create_time_plot(df_plot, "mean_cells", "Mean Cells / Sample", is_log_x = TRUE)
p_markers <- create_time_plot(df_plot, "n_markers", "Number of Markers")
p_samples <- create_time_plot(df_plot, "n_samples", "Number of Samples")
p_pops    <- create_time_plot(df_plot, "n_populations", "Number of Populations")

final_fig_time <- (p_cells + p_markers) / (p_samples + p_pops) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 9),
    plot.tag = element_text(face = "bold", size = 12)
  )

# ------------------------------------------------------------------------------
# 5. SAVE
# ------------------------------------------------------------------------------
ggsave(output_file, plot = final_fig_time, width = 180, height = 160, units = "mm", dpi = 600)

cat(paste0("Success: Figure 4 saved to ", output_file, "\n"))

###
# Usage example:
#Rscript Fig4_plot.r \
#  --perf_input ./ob-blob-metrics/out/performances.tsv \
#  --meta_json ./ob-blob-metrics/out/metric_collectors/metrics_report/dataset_metadata.json \
#  --output ./ob-pipeline-plots/Figure4_Runtime_Metadata.png
###