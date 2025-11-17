## Sonia cebrian camison
## proyecto culex snps
# 11/02/2025
# Script to plot Missing data per individuo of a vcf file

# Load libraries
library(ggplot2)
library(readr)

#  USER VARIABLES
input_file <- "530perex_perex.imiss"         # Path to your input file
main_title <- paste0("Miss gt/indiv ", input_file)  # Main title of the plot
output_file <- paste0("MISS", input_file,".png")   # Name of the output PNG file

data <- read_tsv(input_file)

# Calculate mean F_MISS
mean_fmiss <- mean(data$F_MISS)
sd_fmiss <- sd(data$F_MISS)

# Define upper and lower thresholds
upper_threshold <- mean_fmiss + 2 * sd_fmiss
lower_threshold <- mean_fmiss - 2 * sd_fmiss  # Not always meaningful for missingness, but included


# Create the plot
p <- ggplot(data, aes(x = INDV, y = F_MISS)) +
  geom_point(size = 2) +
  geom_hline(yintercept = mean_fmiss, color = "blue", linetype = "dashed", size = 1) +
  geom_hline(yintercept = upper_threshold, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = lower_threshold, color = "red", linetype = "dashed", size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = main_title,
       subtitle = paste0("Mean = ", round(mean_fmiss, 3),
                         " | Upper Threshold = ", round(upper_threshold, 3),
                         " | Lower Threshold = ", round(lower_threshold, 3)),
       x = "Sample",
       y = "Fraction Missing (F_MISS)")
p
# Save the plot as PNG
ggsave(output_file, plot = p, width = 12, height = 6, dpi = 300)

# Also show plot in RStudio Viewer
print(p)

