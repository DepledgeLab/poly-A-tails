library(tidyverse)
library(readxl)

#setwd("")

# Load the data
Virus <- read_excel("data/250814_HCMV_polya_wTemporalClass.xlsx")

# Rename columns
colnames(Virus) <- c("Sample", "Time", "Median", "Mean", "Mode", "count", "class")

# Fill down Sample and class
Virus <- Virus %>%
  fill(Sample, class, .direction = "down") %>%
  mutate(
    Time = factor(Time, levels = c("24h", "48h", "72h")),
    Time_num = as.numeric(Time)
  )

# Summary statistics per time point
summary_stats <- Virus %>%
  group_by(Time) %>%
  summarise(
    Q1 = quantile(Median, 0.25, na.rm = TRUE),
    Q3 = quantile(Median, 0.75, na.rm = TRUE),
    Median_val = median(Median, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Function to plot a given temporal class
plot_class <- function(class_name, title = class_name, save_file) {
  data_to_plot <- Virus %>%
    filter(class == class_name)
  
  p <- ggplot(data_to_plot, aes(x = Time_num, y = Median, group = Sample)) +
    geom_line(color = "lightyellow4", alpha = 0.6) +
    geom_point(color = "black", size = 1) +
    geom_boxplot(
      aes(group = Time_num),
      fill = "lightyellow",
      alpha = 1,
      width = 0.4,
      position = position_identity()
    ) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("24h", "48h", "72h")) +
    ylim(0, 300) +
    labs(x = "time post infection", y = "median poly(A) length") +
    ggtitle(title) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 25),
      axis.text.x  = element_text(size = 25, color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      axis.title.y = element_text(size = 25),
      axis.text.y  = element_text(size = 25, color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.length = unit(5, "pt"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  print(p)
  ggsave(save_file, plot = p, width = 7.5, height = 5, units = "in")
}

# Example: plot unclassified class ("na")
plot_class(
  class_name = "na",
  title = "unclassified",
  save_file = "HCMV-median-polyA-timecourse-unclassified.pdf"
)

# Example: plot temporal class 7
# plot_class(
#   class_name = "tc_7",
#   title = "temporal class 7",
#   save_file = "HCMV-median-polyA-timecourse-temporal-class-7.pdf"
# )
