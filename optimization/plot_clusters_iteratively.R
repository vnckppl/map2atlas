#!/usr/bin/env Rscript
## * Arguments
args <- commandArgs(trailingOnly = TRUE)
# Location where counts.csv and optimizations.csv live
dir <- args[1]

## * Libraries
library(ggplot2)

## * Input and output files
ifile1 <- file.path(dir, "counts.csv")
ifile2 <- file.path(dir, "optimization.csv")
ofile <- file.path(dir, "optimal_plot.pdf")

## * Load data
data1 <- read.csv(ifile1)
data2 <- read.csv(ifile2)

## * Match data
## ** Create match id for data1
data1$match <- paste(
  as.character(sprintf("%03d", data1$othr1)),
  as.character(sprintf("%03d", data1$othr2)),
  sep = "_"
)
## *** Keep only the cluster count
data1 <- data1[, c("match", "cluster_count")]

## ** Create match id for data2
data2$match <- paste(
  as.character(sprintf("%03d", data2$othr1)),
  as.character(sprintf("%03d", data2$othr2)),
  sep = "_"
)
data <- merge(data2, data1, by = "match", all = TRUE)

## * Prepare data for plotting
data$othr2f <- as.factor(data$othr2)

## ** Variable to highlight combinations of threshold 1 and 2 where the
## proportion outside the mask >= 5%
data$prop_outside_5[data$prop_outside >= 0.05] <-
  data$prop_outside[data$prop_outside >= 0.05]

data_po5 <- data[data$prop_outside >= 0.05, ]

## * Return most optimal combination of thresholds
# Most optimal is where the difference is maximized and the percentage selected
# voxels outside the mask is smaller than 5%.
data_opt <- data[data$prop_outside < 0.05, ]
data_opt <- data_opt[order(data_opt$difference, decreasing = TRUE), ]
thr1_opt <- data_opt$othr1[1]
thr2_opt <- data_opt$othr2[1]
diff_opt <- data$difference[data$othr1 == thr1_opt & data$othr2 == thr2_opt]
my_msg <- paste(
  "Optimal threshold settings: Threshold 1 =", thr1_opt, "; Threshold 2 =",
  thr2_opt
)
print(my_msg)

## * Plot image
p <- ggplot(data, aes(othr1, difference, group = othr2, color = othr2f)) +
  geom_line() +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(size = 14)
  ) +
  xlab("Threshold 1") +
  ylab(
    "Difference in Proportion of Selected Voxels Inside and Outside the Mask"
  ) +
  ggtitle(my_msg) +
  scale_x_continuous(breaks = seq(5, 100, 5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(title = "Threshold 2")) +
  geom_point(
    data = data_po5, aes(othr1, difference),
    shape = 4, size = 2, show.legend = FALSE
  ) +
  annotate(geom = "text", x = thr1_opt, y = diff_opt, label = "o", size = 6)

p

## * Save image
ggsave(ofile, width = 10, height = 10)
