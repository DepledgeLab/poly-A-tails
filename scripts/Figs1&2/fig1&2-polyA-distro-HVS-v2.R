library(tidyverse)
library(ggplot2)
library(biomaRt)

# Mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Paths
#setwd("")

# Load input data (ensure cellular and viral directories are in the same folder as this script)
host_df <- read.table(
  "cellular/callithrix_genome_HVS.polyA.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

virus_df <- read.table(
  "viral/viral_HVS.polyA.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

plot_title <- "HVS, OMK cells, 48 hpi"
plot_file  <- "HVS-callithrix-polyA.pdf"

# Process host
host_clean <- host_df %>%
  filter(qc_tag == "PASS") %>%
  select(contig, polya_length) %>%
  rename(Transcript = contig, PolyA_Length = polya_length) %>%
  mutate(
    Transcript = sub("\\..*", "", Transcript),
    PolyA_Length = as.numeric(PolyA_Length),
    Run = "Callithrix"
  )

# Process virus
virus_clean <- virus_df %>%
  filter(qc_tag == "PASS") %>%
  select(contig, polya_length) %>%
  rename(Transcript = contig, PolyA_Length = polya_length) %>%
  mutate(
    Transcript = sub("\\..*", "", Transcript),
    PolyA_Length = as.numeric(PolyA_Length),
    Run = "HVS"
  )

# Combine
full_polya <- bind_rows(host_clean, virus_clean)

# Mode and summaries (kept for parity with original script)
summary(host_clean)
getmode(round(host_clean$PolyA_Length, 0))

summary(virus_clean)
getmode(round(virus_clean$PolyA_Length, 0))

# Plot
p <- ggplot(full_polya, aes(
  x = Run,
  y = PolyA_Length,
  fill = Run
)) +
  geom_violin(linewidth = 1) +
  geom_boxplot(width = 0.1, linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(
    values = c(Callithrix = "lightblue", HVS = "lightyellow")
  ) +
  ggtitle(plot_title) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, colour = "black"),
    plot.title.position = "plot",
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(linewidth = 0.5, colour = "darkgrey"),
    panel.grid.minor = element_line(linewidth = 0.1, colour = "darkgrey"),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(colour = "black", linewidth = 1),
    axis.line.y = element_line(colour = "black", linewidth = 1),
    axis.title.x = element_text(size = 25, colour = "black"),
    axis.title.y = element_text(size = 25, colour = "black"),
    axis.text.x  = element_text(size = 25, colour = "black"),
    axis.text.y  = element_text(size = 20, colour = "black")
  )

print(p)
ggsave(plot_file, plot = p, width = 5, height = 5, units = "in")
