#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(viridis)
})

# ------------------------------------------------------------------------------
# COMMAND LINE ARGUMENTS
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--conf_input"), type = "character", default = NULL, 
              help = "Path to the per_population_confusion.tsv file", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./", 
              help = "Directory where the individual PNGs will be saved", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$conf_input)) {
  print_help(opt_parser)
  stop("Input file (-i / --conf_input) must be provided.", call. = FALSE)
}

input_file <- opt$conf_input
output_dir <- opt$output_dir

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 1. DATA LOADING & PROCESSING
# ------------------------------------------------------------------------------
df <- read_tsv(input_file, show_col_types = FALSE)

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

df_clean <- df %>%
  mutate(dataset = recode(dataset, !!!name_map)) %>%
  mutate(
    tp = as.numeric(tp),
    fp = as.numeric(fp),
    fn = as.numeric(fn),
    denom = (2 * tp) + fp + fn,
    f1_score = ifelse(denom > 0, (2 * tp) / denom, 0)
  )

# ------------------------------------------------------------------------------
# 2. STYLE DEFINITION
# ------------------------------------------------------------------------------
theme_gb_white <- theme_classic(base_size = 8, base_family = "sans") +
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

# ------------------------------------------------------------------------------
# 3. LOOP & SAVE INDIVIDUALLY
# ------------------------------------------------------------------------------
unique_datasets <- unique(df_clean$dataset)
cat(paste("Found", length(unique_datasets), "datasets. Generating plots...\n"))

for (ds_name in unique_datasets) {
  
  plot_data <- df_clean %>% filter(dataset == ds_name)
  n_pops <- length(unique(plot_data$population_name))
  n_tools <- length(unique(plot_data$model))
  
  p <- ggplot(plot_data, aes(x = population_name, y = f1_score, fill = model)) +
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
    theme_gb_white +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Dynamic Width
  required_width_mm <- (n_pops * n_tools * 5) + 50
  final_width <- min(max(required_width_mm, 100), 400)
  
  safe_name <- str_replace_all(ds_name, "[^A-Za-z0-9]", "_")
  filename <- file.path(output_dir, paste0("Figure2_", safe_name, ".png"))
  
  ggsave(filename, plot = p, width = final_width, height = 100, 
         units = "mm", dpi = 600, device = "png")
  
  cat(paste("Saved:", filename, "\n"))
}

cat("Success: All Figure 2 plots have been generated.\n")

###
# Usage example:
#Rscript Fig2_plot.r \
#  --conf_input ./ob-blob-metrics/out/metric_collectors/metrics_report/per_population_confusion.tsv \
#  --output_dir ./ob-pipeline-plots
###