library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)

# Load------------------
Virus <- read_excel("data/Final-HSV1-polyA-timecourse-wKinetics.xlsx"
)

colnames(Virus) <- c("Sample", "Time", "Median", "Mean", "Mode", "count", "Class")

Virus <- Virus %>%
  fill(Sample, Class, .direction = "down")

# Subsets---------------
Virus_sel <- Virus %>%
  filter(Sample %in% c("RL2-1","RS1-1","UL1-1","UL46-1","UL47-1",
                       "UL15-1","UL24-1","UL52-1","US1-1"))

# Function for summary per timepoint --------------------------------------
summ_stats <- function(df, t) {
  df %>%
    filter(Time == t) %>%
    summarise(
      Q1 = quantile(Median, 0.25, na.rm = TRUE),
      Q3 = quantile(Median, 0.75, na.rm = TRUE),
      Median_val = median(Median, na.rm = TRUE)
    ) %>%
    mutate(Time = t)
}

summary_stats <- bind_rows(
  summ_stats(Virus_sel, "3h"),
  summ_stats(Virus_sel, "6h"),
  summ_stats(Virus_sel, "12h")
)

# Factor ordering-------
time_levels <- c("3h", "6h", "12h")

Virus <- Virus %>%
  mutate(Time_num = as.numeric(factor(Time, levels = time_levels)))

summary_stats <- summary_stats %>%
  mutate(Time = factor(Time, levels = time_levels),
         Time_num = as.numeric(Time))

# Split by kinetic class
Virus_classes <- Virus %>%
  group_split(Class, .keep = TRUE) %>%
  set_names(unique(Virus$Class))

Virus_alpha  <- Virus_classes[["alpha"]]
Virus_beta   <- Virus_classes[["beta"]]
Virus_gamma  <- Virus_classes[["gamma"]]
Virus_gamma1 <- Virus_classes[["gamma1"]]
Virus_gamma2 <- Virus_classes[["gamma2"]]

# Choose which class to plot ----------------------------------------------
df_plot <- Virus_gamma2
save_path <- "HSV1-median-polyA-timecourse-gamma2.pdf"

# Plot------------------
p <- ggplot(df_plot, aes(x = Time_num, y = Median, group = Sample)) +
  geom_line(color = "lightyellow4", alpha = 0.6) +
  geom_point(color = "black", size = 1) +
  geom_boxplot(
    aes(group = Time_num),
    fill = "lightyellow",
    width = 0.4,
    position = position_identity()
  ) +
  scale_x_continuous(breaks = 1:3, labels = c("3", "6", "12")) +
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
ggsave(save_path, plot = p, width = 7.5, height = 5, units = "in")
