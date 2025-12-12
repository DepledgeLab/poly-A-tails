library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)

# Load data ---------------------------------------------------------------
Virus <- read_excel("data/Final-HSV1-polyA-timecourse.xlsx")

colnames(Virus) <- c("Sample", "Time", "Median", "Mean", "Mode", "count")

# Fill down sample names if they are in the first row only -----------------
Virus <- Virus %>%
  fill(Sample, .direction = "down")

# Convenience sets ---------------------------------------------------------
Virus_sel <- Virus %>%
  filter(Sample %in% c("RL2-1", "RS1-1", "US12-1"))

Virus_all <- Virus %>%
  filter(!Sample %in% "RL2-1")

# Summaries for 3h / 6h / 12h ----------------------------------------------
summ_stats <- function(df, time_point) {
  df %>%
    filter(Time == time_point) %>%
    summarise(
      Q1 = quantile(Median, 0.25, na.rm = TRUE),
      Q3 = quantile(Median, 0.75, na.rm = TRUE),
      Median_val = median(Median, na.rm = TRUE)
    ) %>%
    mutate(Time = time_point)
}

summary_stats <- bind_rows(
  summ_stats(Virus_sel, "3h"),
  summ_stats(Virus_sel, "6h"),
  summ_stats(Virus_sel, "12h")
)

# Factor time in the correct order -----------------------------------------
time_levels <- c("3h", "6h", "12h")

Virus_sel <- Virus_sel %>%
  mutate(Time_num = as.numeric(factor(Time, levels = time_levels)))

Virus <- Virus %>%
  mutate(Time_num = as.numeric(factor(Time, levels = time_levels)))

summary_stats <- summary_stats %>%
  mutate(Time = factor(Time, levels = time_levels),
         Time_num = as.numeric(Time))

# Plot 1: all viral mRNAs --------------------------------------------------
p <- ggplot(Virus, aes(x = Time_num, y = Median, group = Sample)) +
  geom_line(color = "lightyellow4", alpha = 0.6) +
  geom_point(color = "black", size = 1) +
  geom_boxplot(
    aes(group = Time_num),
    fill = "lightyellow",
    width = 0.4,
    position = position_identity()
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("3", "6", "12")
  ) +
  ylim(0, 200) +
  labs(
    x = "hours post infection (hpi)",
    y = "median poly(A) tail length "
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 25),
    axis.text.x  = element_text(size = 25, color = "black"),
    axis.title.y = element_text(size = 25),
    axis.text.y  = element_text(size = 25, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

print(p)

ggsave("HSV1-median-polyA-timecourse.pdf",
  plot = p, width = 7.5, height = 5
)

# Plot 2: selected viral mRNAs ---------------------------------------------
q <- ggplot(Virus_sel, aes(x = Time_num, y = Median, group = Sample)) +
  geom_line(color = "goldenrod4", alpha = 0.6, size = 1.5) +
  geom_point(color = "black", size = 1) +
  geom_boxplot(
    aes(group = Time_num),
    fill = "goldenrod1",
    width = 0.4,
    position = position_identity()
  ) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("3", "6", "12")
  ) +
  ylim(0, 200) +
  labs(
    x = "hours post infection (hpi)",
    y = "median poly(A) tail length "
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 25),
    axis.text.x  = element_text(size = 25, color = "black"),
    axis.title.y = element_text(size = 25),
    axis.text.y  = element_text(size = 25, color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

print(q)

ggsave("HSV1-median-polyA-timecourse-select-new.pdf",
  plot = q, width = 7.5, height = 5
)
