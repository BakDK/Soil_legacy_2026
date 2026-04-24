# -------------------------------------------------------------------------
# # Plotting figure 5. 
# The Random forest modelling was performed in another script and output files stored in: "outputs/legacy/block/models_and_shap"
# -------------------------------------------------------------------------
# 
library(tidyverse)
library(patchwork)
library(here) # Crucial for robust file paths

# 1. Define Paths relative to the Git Root
# 'here()' automatically finds the 'Future_cropping' folder
input_dir  <- here("outputs", "legacy", "legacy_non_rarefied_block_train13_test24", "02_models_and_shap")
output_dir <- here("Plots")

# Ensure output directory exists
if(!dir.exists(output_dir)) dir.create(output_dir)

# -------------------------------------------------------------------------
# DATA LOADING
# -------------------------------------------------------------------------

# Define DAS values
das_values <- c(11, 14, 17, 31, 158, 187, 215, 227, 256, 318)

# Construct full file paths using the robust input_dir
file_names <- sprintf("rf_performance_metrics_test_%dDAS.csv", das_values)
full_file_paths <- file.path(input_dir, file_names)

# Read and combine data
full_data <- map_dfr(full_file_paths, function(file) {
  
  # Check if file exists to prevent pipeline crashes
  if (!file.exists(file)) {
    warning(paste("File not found:", file))
    return(NULL)
  }
  
  # Extract DAS number
  das_val <- as.numeric(str_extract(basename(file), "\\d+(?=DAS)"))
  
  # Read CSV with explicit column names
  read_csv(file, col_names = c("metric", "Mean", "Std_Dev"), 
           skip = 1, show_col_types = FALSE) %>%
    mutate(DAS_Num = das_val)
})

full_data$sem<-full_data$Std_Dev/sqrt(20) #20 RF runs

# -------------------------------------------------------------------------
# DATA PREPARATION
# -------------------------------------------------------------------------

# Prepare Panel a (AUC)
df_auc <- full_data %>% 
  filter(metric %in% c("auc_M1P1", "auc_N1K1", "auc_N1P2K2")) %>%
  mutate(
    Day_Label = paste0(DAS_Num, "DAS"),
    Day_Factor = factor(Day_Label, levels = paste0(das_values, "DAS")),
    Fertility = metric
  )

# Prepare Panel b (Kappa)
df_kappa <- full_data %>%
  filter(metric == "kappa") %>%
  mutate(
    Day_Factor = factor(DAS_Num, levels = das_values)
  )

# -------------------------------------------------------------------------
# PLOTTING
# -------------------------------------------------------------------------

# Panel a: AUC Lines
plot_a <- ggplot(df_auc, aes(x = Day_Factor, y = Mean, colour = Fertility, group = Fertility)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - Std_Dev, ymax = Mean + Std_Dev), width = 0.2) +
  geom_line() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = c(6.5, 7.5), linetype = "dashed", linewidth = 0.7, color = "Black") +
  scale_y_continuous(name = "AUC", breaks = c(0.5, 0.6, 0.8, 1.0)) +
  scale_color_manual(
    values = c("#E69F00","#999999" , "#56B4E9"), 
    labels = c("M1P1", "N1K1", "N1P2K2")
  ) +
  scale_x_discrete(labels = as.character(das_values)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.position = "inside",
    legend.position.inside = c(0.24, 0.1), 
    legend.direction = "horizontal",
    axis.title.x = element_blank() 
  ) +
  labs(colour = "Fertility Level")

# Panel b: Kappa Bars
plot_b <- ggplot(df_kappa, aes(x = Day_Factor, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean - Std_Dev, ymax = Mean + Std_Dev), width = 0.25, linewidth = 0.8) +
  labs(x = "Days After Sowing (DAS)", y = "Kappa") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.line = element_line(color = "black")
  )

# Combine Plots
combined_plot <- (plot_a / plot_b) + 
  plot_annotation(tag_levels = 'A') &  # This will label the panels as "a" and "b"
  theme(plot.tag = element_text(size = 20, face = "bold"))

combined_plot
# -------------------------------------------------------------------------
# SAVING
# -------------------------------------------------------------------------

# Save to the specific 'Figures' folder defined at the top
output_file <- file.path(output_dir, "Figure5_nonrare_13.png")

ggsave(output_file, plot = combined_plot, width = 10, height = 12, dpi = 300)

message(paste("Plot saved successfully to:", output_file))
