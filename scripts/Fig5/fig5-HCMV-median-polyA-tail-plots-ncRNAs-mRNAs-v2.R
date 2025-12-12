library(tidyverse)
library(readxl)

#setwd("")


#Load and prepare data
Virus <- read_excel("data/250814_HCMV_polya.xlsx")
#Virus <- read_excel("data/250814_HCMV_ncRNA.xlsx")

colnames(Virus) <- c("Sample", "Time", "Median", "Mean", "Mode", "count")

Virus <- Virus %>%
  fill(Sample, .direction = "down") %>%
  mutate(
    Time = factor(Time, levels = c("24h", "48h", "72h")),
    Time_num = as.numeric(Time)
  )


#Summary statistics per time point
summary_stats <- Virus %>%
  group_by(Time) %>%
  summarise(
    Q1 = quantile(Median, 0.25, na.rm = TRUE),
    Q3 = quantile(Median, 0.75, na.rm = TRUE),
    Median_val = median(Median, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)


#Plot
p <- ggplot(Virus, aes(x = Time_num, y = Median, group = Sample)) +
  geom_line(color = "lightyellow4", alpha = 0.6) +
  geom_point(color = "black", size = 1) +
  geom_boxplot(
    aes(group = Time_num),
    fill = "lightyellow",
    alpha = 1,
    width = 0.4,
    position = position_identity()
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = c("24h", "48h", "72h")
  ) +
  ylim(0, 300) +
  labs(
    x = "time post infection",
    y = "median poly(A) length "
  ) +
  # ggtitle("HCMV ncRNAs") +
  ggtitle("HCMV mRNAs") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 25),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    
    axis.title.y = element_text(size = 25),
    axis.text.y = element_text(size = 25, color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    
    axis.ticks.length = unit(5, "pt"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

print(p)

save <- "HCMV-median-polyA-timecourse-mRNA.pdf"
#save <- "HCMV-median-polyA-timecourse-ncRNA.pdf"

ggsave(save, plot = p, width = 7.5, height = 5, units = "in")
