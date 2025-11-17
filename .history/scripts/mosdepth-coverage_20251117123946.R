#!/usr/bin/env Rscript

# =============================================================================
# Coverage analysis and plotting from mosdepth outputs (current directory version)
# Usage: Rscript mosdepth-coverage.R
# =============================================================================

library(tidyverse)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Auto-detect sample names from global.dist filenames
# -----------------------------------------------------------------------------
file_list <- list.files(pattern = "\\.merged\\.mosdepth\\.global\\.dist\\.txt$")

if (length(file_list) == 0) {
  stop("No mosdepth global distribution files found.")
}

# Extract sample names from filenames
samples <- gsub("\\.merged\\.mosdepth\\.global\\.dist\\.txt$", "", file_list)

# -----------------------------------------------------------------------------
# 2. Load GLOBAL coverage data
# -----------------------------------------------------------------------------
global_data_list <- lapply(samples, function(s) {
  file_gdist <- paste0(s, ".merged.mosdepth.global.dist.txt")

  if (!file.exists(file_gdist)) {
    warning(paste("Missing file:", file_gdist))
    return(NULL)
  }

  df <- read_tsv(file_gdist,
                 col_names = c("chrom", "coverage", "fraction"),
                 comment = "#") %>%
    mutate(sample = s)

  return(df)
})

global_data_list <- global_data_list[!sapply(global_data_list, is.null)]
global_data <- bind_rows(global_data_list)

# -----------------------------------------------------------------------------
# 3. Load REGIONS coverage data
# -----------------------------------------------------------------------------
regions_data_list <- lapply(samples, function(s) {
  file_regions <- paste0(s, ".merged.regions.bed.gz")

  if (!file.exists(file_regions)) {
    warning(paste("Missing regions file for:", s))
    return(NULL)
  }

  df <- read_tsv(file_regions,
                 col_names = c("chrom", "start", "end", "mean_coverage"),
                 comment = "#") %>%
    mutate(sample = s,
           mid = (start + end) / 2)

  return(df)
})

regions_data_list <- regions_data_list[!sapply(regions_data_list, is.null)]
regions_data <- bind_rows(regions_data_list)

# -----------------------------------------------------------------------------
# 4. PLOT 1: Global Coverage Distribution
# -----------------------------------------------------------------------------
p <- "coverage_distribution_all_samples.png"

dist <- ggplot(global_data, aes(x = coverage, y = fraction, color = sample)) +
  geom_line(linewidth = 1, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Global Coverage Distribution",
    x = "Coverage",
    y = "Proportion of Genome ≥ Coverage"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(limits = c(0, 140), breaks = seq(0, 140, 10)) +
  scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1.00, 0.10))

ggsave(p, device = "png", units = "in", plot = dist, width = 6, height = 4)

# -----------------------------------------------------------------------------
# 5. PLOT 2: Per-Chromosome Coverage
# -----------------------------------------------------------------------------
chroms <- unique(regions_data$chrom)

chrom_plots <- lapply(chroms, function(chr) {
  chr_data <- regions_data %>% filter(chrom == chr)

  ggplot(chr_data, aes(x = mid / 1e6, y = mean_coverage, color = sample)) +
    geom_line(linewidth = 0.8, alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste("Coverage Along Chromosome", chr),
      x = "Position (Mbp)",
      y = "Mean Depth"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(0, 50))
})

combined_plot <- wrap_plots(chrom_plots, ncol = 1)

f_combined <- "coverage_chromosomes_all_samples.png"
ggsave(f_combined, plot = combined_plot, width = 8, height = 4 * length(chroms), units = "in")