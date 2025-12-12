library(tidyverse)
library(ggplot2)


# Helper function for repeated preprocessing
clean_polya <- function(df, contig_pattern, run_label = NULL) {
  df %>%
    filter(str_detect(contig, contig_pattern)) %>%
    select(contig, polya_length) %>%
    rename(Transcript = contig, PolyA_Length = polya_length) %>%
    mutate(
      Transcript = str_replace(Transcript, "\\..*", ""),
      PolyA_Length = as.numeric(PolyA_Length),
      Run = run_label
    )
}


# Load TB40 viral dataset
CTRL1 <- read.table(
  "fig3_data/TB40_72h.TB40txome.polyA.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
) %>%
  filter(qc_tag == "PASS")

# Process viral transcripts of interest
RNA4_9L <- clean_polya(CTRL1, "^mRNA.RNA4.9L",   "ncRNA 4.9")
RNA1_2  <- clean_polya(CTRL1, "^RNA1.2-1",        "ncRNA 1.2")
RNA2_7  <- clean_polya(CTRL1, "^RNA2.7-1",        "ncRNA 2.7")
RNA5_0  <- clean_polya(CTRL1, "^RNA5.0-1",        "RNA 5.0 precursor")

# List of contigs to exclude for "remaining" group
exclude_contigs <- c(
  "mRNA.RNA4.9L::TB40E-cloneBAC4:126772-131672",
  "RNA1.2-1::TB40E-cloneBAC4:39120-40144",
  "RNA2.7-1::TB40E-cloneBAC4:35218-38565",
  "RNA5.0-1::TB40E-cloneBAC4:188496-194174"
)

remain <- CTRL1 %>%
  filter(!contig %in% exclude_contigs) %>%
  select(contig, polya_length) %>%
  rename(Transcript = contig, PolyA_Length = polya_length) %>%
  mutate(
    Transcript = str_replace(Transcript, "\\..*", ""),
    Run = "HCMV mRNAs"
  )

# Load cellular ncRNAs (MALAT1 and NEAT1)
cell_ncRNAs <- read.table(
  "fig3_data/TB40_72h.GRCh38.tx.polyA.tsv",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)

colnames(cell_ncRNAs) <- c(
  "readname","contig","position","leader_start","adapter_start",
  "polya_start","transcript_start","read_rate","polya_length","qc_tag"
)

cell_ncRNAs <- cell_ncRNAs %>%
  filter(qc_tag == "PASS")

MALAT <- clean_polya(cell_ncRNAs, "MALAT", "MALAT")
NEAT  <- clean_polya(cell_ncRNAs, "NEAT1", "NEAT")

ncRNA <- bind_rows(MALAT, NEAT)


# Combine everything for plotting
full_polya <- bind_rows(
  remain,
  RNA4_9L,
  RNA1_2,
  RNA2_7,
  NEAT       # excluding RNA5.0 precursor as before
) %>%
  mutate(
    Transcript = factor(Transcript),
    PolyA_Length = as.numeric(PolyA_Length)
  )


# Calculate modes (same logic, cleaner code)
getmode <- function(v) {
  v <- round(as.numeric(v), 0)
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(NEAT$PolyA_Length)
getmode(RNA1_2$PolyA_Length)
getmode(remain$PolyA_Length)
getmode(RNA2_7$PolyA_Length)
getmode(RNA4_9L$PolyA_Length)

# Plot
title <- "HCMV, NHDF, MOI=3, 72hpi"

p <- ggplot(full_polya, aes(
  x = factor(Run, levels = c("HCMV mRNAs","ncRNA 4.9","ncRNA 1.2","ncRNA 2.7","NEAT"),
             labels = c("HCMV mRNAs","ncRNA4.9","ncRNA1.2","ncRNA2.7","NEAT1")),
  y = PolyA_Length,
  fill = Run
)) +
  geom_violin(linewidth = 1) +
  geom_boxplot(width = 0.1, linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(values = c(
    "HCMV mRNAs"="lightyellow", "ncRNA 4.9"="coral1",
    "ncRNA 1.2"="coral1", "ncRNA 2.7"="coral1", "NEAT"="lightblue1"
  )) +
  ggtitle(title) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, colour = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(linewidth = 0.5, colour = "darkgrey"),
    panel.grid.minor = element_line(linewidth = 0.1, colour = "darkgrey"),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    axis.title = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 20)
  )

print(p)

ggsave(
  "HCMV-ncRNA-analysis.pdf",
  plot = p, width = 12.5, height = 5, units = "in"
)
