# =============================================================
# 02_clonotype_enrichment.R
# Clonotype Sharing & Enrichment Analysis — TCRβ in T1D vs CTRL
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

# ---- LOAD METADATA ------------------------------------------
metadata <- read_csv("metadata.csv")

meta_tc <- metadata %>%
  filter(study_group_description %in% c("T1D", "CTRL")) %>%
  select(filename, study_group_description)

# ---- EXTRACT UNIQUE CLONOTYPES PER SAMPLE -------------------
# For each individual, extract unique clonotypes defined by:
# junction_aa (CDR3 amino acid sequence) + v_call + j_call

presence_tbl <- meta_tc %>%
  mutate(
    clonotypes = map(
      filename,
      ~ read_tsv(
          file.path(folder_path, .x),
          col_types = cols_only(
            junction_aa = col_character(),
            v_call      = col_character(),
            j_call      = col_character()
          ),
          progress = FALSE
        ) %>%
        distinct(junction_aa, v_call, j_call) %>%
        transmute(clonotype = paste(junction_aa, v_call, j_call, sep = "|")) %>%
        pull(clonotype)
    )
  ) %>%
  unnest(clonotypes) %>%
  rename(clonotype = clonotypes)

# ---- BUILD PRESENCE TABLE -----------------------------------
# Count how many individuals in each group carry each clonotype
# (presence/absence per individual, not abundance)

clono_counts <- presence_tbl %>%
  count(clonotype, study_group_description, name = "n_present") %>%
  pivot_wider(
    names_from  = study_group_description,
    values_from = n_present,
    values_fill = 0
  )

# ---- SHARED CLONOTYPES --------------------------------------
# Clonotypes present in at least one T1D AND one CTRL individual

shared_clonotypes <- clono_counts %>%
  filter(T1D >= 1, CTRL >= 1)

cat("Total shared clonotypes (T1D ∩ CTRL):", nrow(shared_clonotypes), "\n")

# ---- ENRICHMENT CALCULATION ---------------------------------
# Compute prevalence per group and log2 enrichment ratio
# log2 > 0 = enriched in T1D
# log2 < 0 = enriched in CTRL

group_sizes <- meta_tc %>%
  count(study_group_description)

N_CTRL <- group_sizes$n[group_sizes$study_group_description == "CTRL"]
N_T1D  <- group_sizes$n[group_sizes$study_group_description == "T1D"]

cat("Group sizes — T1D:", N_T1D, "| CTRL:", N_CTRL, "\n")

enriched_clonotypes <- shared_clonotypes %>%
  mutate(
    prev_T1D        = T1D  / N_T1D,
    prev_CTRL       = CTRL / N_CTRL,
    enrichment_log2 = log2(prev_T1D / prev_CTRL)
  )

# ---- SUMMARY ------------------------------------------------
enrichment_summary <- enriched_clonotypes %>%
  summarise(
    n_total             = n(),
    n_enriched_T1D      = sum(enrichment_log2 > 0),
    n_enriched_CTRL     = sum(enrichment_log2 < 0),
    median_enrichment   = median(enrichment_log2),
    mean_enrichment     = mean(enrichment_log2)
  )

print(enrichment_summary)

# ---- SAVE OUTPUT --------------------------------------------
write_csv(enriched_clonotypes, "results/enriched_clonotypes.csv")

# ---- VISUALIZATION ------------------------------------------
# Distribution of log2 enrichment scores

p_enrichment <- ggplot(enriched_clonotypes,
                       aes(x = enrichment_log2)) +
  geom_histogram(binwidth = 0.1, fill = "#1B9E77", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title    = "Distribution of Clonotype Enrichment Scores",
    subtitle = "log2(prevalence T1D / prevalence CTRL) for shared clonotypes",
    x        = "log2 Enrichment Score",
    y        = "Number of Clonotypes"
  ) +
  theme_minimal(base_size = 12)

ggsave("results/clonotype_enrichment_distribution.png",
       p_enrichment, width = 10, height = 6, dpi = 300)

# Top clonotypes enriched in T1D
top_T1D <- enriched_clonotypes %>%
  arrange(desc(enrichment_log2)) %>%
  head(20)

# Top clonotypes enriched in CTRL
top_CTRL <- enriched_clonotypes %>%
  arrange(enrichment_log2) %>%
  head(20)

write_csv(top_T1D,  "results/top20_enriched_T1D.csv")
write_csv(top_CTRL, "results/top20_enriched_CTRL.csv")

cat("Done. Results saved to results/\n")
