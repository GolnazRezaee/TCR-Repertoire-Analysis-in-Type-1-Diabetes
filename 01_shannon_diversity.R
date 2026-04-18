# =============================================================
# 01_shannon_diversity.R
# Shannon Diversity Analysis — TCRβ Repertoire in T1D vs CTRL
# DTU Erasmus Internship 2025-2026
# Supervisor: Prof. Leon Eyrich Jessen
# =============================================================
# NOTE: Raw data not included in this repository.
# Data source: Public AIRR-seq dataset (908 TCRβ repertoire files)
# Place TSV files in data/ and metadata in metadata.csv to reproduce.
# =============================================================

library(tidyverse)

# ---- SETTINGS -----------------------------------------------
folder_path <- "data/"   # replace with your local path to TSV files

# ---- FUNCTION -----------------------------------------------
shannon_diversity_norm <- function(x, w) {
  x_expanded <- rep(x, w)
  x_counts   <- table(x_expanded)
  p_x        <- x_counts / sum(x_counts)
  n          <- length(p_x)
  H          <- -sum(p_x * log2(p_x))
  H_norm     <- H / log2(n)
  return(H_norm)
}

# ---- LOAD FILES ---------------------------------------------
files <- list.files(
  path      = folder_path,
  pattern   = "\\.tsv$",
  full.names = TRUE
)

# ---- COMPUTE DIVERSITY --------------------------------------
results <- data.frame(
  file_name          = character(),
  shannon_diversity  = numeric(),
  n_templates        = integer(),
  n_unique_sequences = integer(),
  stringsAsFactors   = FALSE
)

for (i in seq_along(files)) {
  data <- read_tsv(files[i], show_col_types = FALSE)
  x    <- data$junction_aa
  w    <- data$templates
  H    <- shannon_diversity_norm(x, w)

  results[i, ] <- list(
    file_name          = basename(files[i]),
    shannon_diversity  = H,
    n_templates        = sum(w),
    n_unique_sequences = length(unique(x))
  )

  if (i %% 50 == 0) cat("Processed", i, "out of", length(files), "\n")
}

write_csv(results, "results/shannon_diversity.csv")

# ---- MERGE METADATA -----------------------------------------
metadata    <- read_csv("metadata.csv")
merged_data <- results %>%
  left_join(metadata, by = c("file_name" = "filename"))

data_filtered <- merged_data %>%
  filter(
    study_group_description %in% c("T1D", "CTRL"),
    !is.na(age),
    !is.na(shannon_diversity)
  )

# ---- SUMMARY STATISTICS -------------------------------------
summary_stats <- data_filtered %>%
  group_by(study_group_description) %>%
  summarise(
    n                = n(),
    mean_diversity   = mean(shannon_diversity, na.rm = TRUE),
    sd_diversity     = sd(shannon_diversity,   na.rm = TRUE),
    median_diversity = median(shannon_diversity, na.rm = TRUE),
    q25              = quantile(shannon_diversity, 0.25, na.rm = TRUE),
    q75              = quantile(shannon_diversity, 0.75, na.rm = TRUE),
    .groups          = "drop"
  )
print(summary_stats)

# ---- WILCOXON TEST ------------------------------------------
wilcox_result <- wilcox.test(
  shannon_diversity ~ study_group_description,
  data  = data_filtered,
  exact = FALSE
)
cat("Wilcoxon p-value:", wilcox_result$p.value, "\n")

# ---- VISUALIZATIONS -----------------------------------------
colors <- c("T1D" = "#D95F02", "CTRL" = "#1B9E77")

# Shannon vs Age
p_age <- ggplot(data_filtered,
                aes(x = age, y = shannon_diversity,
                    color = study_group_description)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", span = 1.0, se = TRUE,
              alpha = 0.2, linewidth = 1.2) +
  scale_color_manual(values = colors) +
  labs(title = "Shannon Diversity vs Age",
       x = "Age (years)", y = "Normalized Shannon Diversity",
       color = "Group") +
  theme_minimal(base_size = 12)

ggsave("results/shannon_vs_age.png", p_age, width = 10, height = 6, dpi = 300)

# Shannon vs Age by Sex
data_filtered_sex <- data_filtered %>%
  filter(sex %in% c("F", "M"))

p_sex <- ggplot(data_filtered_sex,
                aes(x = age, y = shannon_diversity,
                    color = study_group_description)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", span = 1.0, se = TRUE,
              alpha = 0.2, linewidth = 1.2) +
  scale_color_manual(values = colors) +
  facet_wrap(~ sex,
             labeller = labeller(sex = c("F" = "Female", "M" = "Male"))) +
  labs(title = "Shannon Diversity vs Age by Sex",
       x = "Age (years)", y = "Normalized Shannon Diversity",
       color = "Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("results/shannon_vs_age_by_sex.png", p_sex, width = 12, height = 6, dpi = 300)
