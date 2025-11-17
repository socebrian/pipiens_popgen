#!/usr/bin/env Rscript
# depth_plot.R  – make depth histograms for a single chromosome
# Usage:  Rscript depth_plot.R <chrom> <ldepth file> <ldepth.mean file> [prefix]

## ---- 1. read the command–line arguments ---------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("Usage: Rscript depth_plot.R <chrom> <ldepth file> <ldepth.mean file> [prefix]")
}

chrom              <- as.integer(args[1])
ldepth_file        <- args[2]
ldepth_mean_file   <- args[3]
prefix             <- ifelse(length(args) >= 4, args[4], tools::file_path_sans_ext(basename(ldepth_file)))

## ---- 2. build the four output names -------------------------------------
name1 <- sprintf("%s_%s_dp_total.png",     prefix, chrom)
name2 <- sprintf("%s_%s_dp_cut.png",       prefix, chrom)
name3 <- sprintf("%s_%s_meandp_total.png", prefix, chrom)
name4 <- sprintf("%s_%s_meandp_cut.png",   prefix, chrom)

## ---- 3. read the tables --------------------------------------------------
site_depth       <- read.table(ldepth_file,       header = TRUE, sep = "\t")
site_mean_depth  <- read.table(ldepth_mean_file,  header = TRUE, sep = "\t")

## ---- 4. plotting code  ------------------------------------
# Check the  structure of the data
str(site_depth)
str(site_mean_depth)

#creat values and explore distribution
mean_depth <- mean(site_depth$SUM_DEPTH)
double_mean_depth <- 2 * mean_depth
max_c1 <- max(site_depth$SUM_DEPTH)
median <- median(site_depth$SUM_DEPTH)

# Plot the histogram, total values
png(name1, width = 800, height = 600)  # Open PNG device with specified size
depthsite_total <- hist(site_depth$SUM_DEPTH,
     breaks = 50,                   # Number of bins
     col = "skyblue",               # Color of the bars
     main = paste0("Depth per site, chromosome", chrom),
     xlab = "SUM_DEPTH",
     ylab = "Frequency",
     border = "black")
abline(v = mean_depth, col = "red", lwd = 2, lty = 2) # Dashed red line
abline(v = median, col = "black", lwd = 2, lty = 2) # Dashed red line
legend("topright",                                       # Position of the legend
       legend = c("Mean", "Median"),                    # Labels for the lines
       col = c("red", "black"),                         # Colors of the lines
       lty = 2,                                         # Line type (dashed)
       lwd = 2,                                         # Line width
       bty = "n")
depthsite_total
dev.off()



#plot histogram cutting on 2*mean
### ESTA ES LA QUE ENVIO DE REPORT
depthsite_dmean <- site_depth[site_depth$SUM_DEPTH < double_mean_depth, ]
png(name2, width = 800, height = 600)  # Open PNG device with specified size
depthsite_dmean_plot <- hist(depthsite_dmean$SUM_DEPTH,
                            breaks = 80,                   # Number of bins
                            col = "skyblue",               # Color of the bars
                            main = paste0("chromosome", chrom),
                            xlab = "Total depth per site (cut at 2*mean)",
                            ylab = "Nº sites",
                            border = "black",
                            xaxp = c(0, 750, 15))
abline(v = mean_depth, col = "red", lwd = 2, lty = 2) # Dashed red line
abline(v = median, col = "black", lwd = 2, lty = 2) # Dashed red line
legend("topright",                                       # Position of the legend
       legend = c("Mean", "Median"),                    # Labels for the lines
       col = c("red", "black"),                         # Colors of the lines
       lty = 2,                                         # Line type (dashed)
       lwd = 2,                                         # Line width
       bty = "n")                                       # No box around the legend


dev.off()


###plot of mean depth, chrom 1 ----

#creat values and explore distribution
mean_depth2 <- mean(site_mean_depth$MEAN_DEPTH)
double_mean_depth2 <- 2 * mean_depth2
max <- max(site_mean_depth$MEAN_DEPTH)
median <- median(site_mean_depth$MEAN_DEPTH)

# Plot the histogram, total values
png(name3, width = 800, height = 600)  
depthmean_total <- hist(site_mean_depth$MEAN_DEPTH,
                            breaks = 50,                   # Number of bins
                            col = "skyblue",               # Color of the bars
                            main = paste0("Mean depth per site, chromosome", chrom),
                            xlab = "SUM_DEPTH",
                            ylab = "Frequency",
                            border = "black")
abline(v = mean_depth2, col = "red", lwd = 2, lty = 2) # Dashed red line
abline(v = median, col = "black", lwd = 2, lty = 2) # Dashed red line
legend("topright",                                       # Position of the legend
       legend = c("Mean", "Median"),                    # Labels for the lines
       col = c("red", "black"),                         # Colors of the lines
       lty = 2,                                         # Line type (dashed)
       lwd = 2,                                         # Line width
       bty = "n")
depthmean_total
dev.off()

### ESTA ES LA QUE ENVIO DE REPORT
site_mean_depth_dmean <- site_mean_depth[site_mean_depth$MEAN_DEPTH < double_mean_depth2, ]
png(name4, width = 800, height = 600)  # Open PNG device with specified size
depthsite_dmean <- hist(site_mean_depth_dmean$MEAN_DEPTH,
                            breaks = 16,                   # Number of bins
                            col = "skyblue",               # Color of the bars
                            main = paste0("chromosome ", chrom),
                            xlab = "Mean depth per site (cut at 2*mean)",
                            ylab = "Nº sites",
                            border = "black",
                            xaxp = c(0, 8, 8))
abline(v = mean_depth2, col = "red", lwd = 2, lty = 2) # Dashed red line
abline(v = median, col = "black", lwd = 2, lty = 2) # Dashed red line
legend("topright",                                       # Position of the legend
       legend = c("Mean", "Median"),                    # Labels for the lines
       col = c("red", "black"),                         # Colors of the lines
       lty = 2,                                         # Line type (dashed)
       lwd = 2,                                         # Line width
       bty = "n")                                       # No box around the legend


dev.off()