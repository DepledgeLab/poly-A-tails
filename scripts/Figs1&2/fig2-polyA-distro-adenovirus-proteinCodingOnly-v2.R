library(tidyverse)
library(ggplot2)
library(biomaRt)
library(stringr)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#setwd("")


# Edit as needed
plot_title <- "hAdv41, A549, MOI=10"
plot_file  <- "Adenovirus-timecourse-polyA-protCodingonly.pdf"


# Virus processing function
process_virus <- function(path, label) {
  read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(qc_tag == "PASS") %>%
    select(contig, polya_length) %>%
    mutate(
      contig = sub("\\..*", "", contig),
      polya_length = as.numeric(polya_length),
      Run = label
    )
}

# Load input data (ensure viral directory is in the same folder as this script)
t12 <- process_virus("viral/viral_AdV41-A549-12h-1.Ad41.polyA.tsv", "t12")
t24 <- process_virus("viral/viral_AdV41-A549-24h-1.Ad41.polyA.tsv", "t24")
t48 <- process_virus("viral/viral_AdV41-A549-48h-1.Ad41.polyA.tsv", "t48")


# Human processing function
process_human <- function(path, label, mart) {
  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(contig_no_version = sub("\\..*", "", contig))
  
  enst_ids <- unique(df$contig_no_version)
  
  annot <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                   "external_gene_name", "transcript_biotype"),
    filters = "ensembl_transcript_id",
    values = enst_ids,
    mart = mart
  )
  
  df %>%
    left_join(annot, by = c("contig_no_version" = "ensembl_transcript_id")) %>%
    filter(qc_tag == "PASS") %>%
    filter(!grepl("Mt_\\|?$", transcript_biotype)) %>%
    filter(
      !is.na(transcript_biotype),
      !grepl(
        "polymorphic_pseudogene|artifact|transcribed_unitary_pseudogene|scaRNA|snRNA|Mt_tRNA|ribozyme|misc_RNA|rRNA|snoRNA|non_stop_decay|processed_transcript|scRNA|protein_coding_LoF|processed_pseudogene|unprocessed_pseudogene|nonsense_mediated_decay|retained_intron|lncRNA",
        transcript_biotype
      ),
      !grepl("^MT-", external_gene_name, ignore.case = TRUE)
    ) %>%
    select(contig = contig_no_version, polya_length) %>%
    mutate(
      contig = sub("\\..*", "", contig),
      polya_length = as.numeric(polya_length),
      Run = label
    )
}

# Single Ensembl connection (no repeated connections needed)
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 104
)

# Load input data (ensure cellular directory is in the same folder as this script)
h12 <- process_human("cellular/human_txome_AdV41-A549-12h-1.GRCh38.polyA.tsv", "h12", ensembl)
h24 <- process_human("cellular/human_txome_AdV41-A549-24h-1.GRCh38.polyA.tsv", "h24", ensembl)
h48 <- process_human("cellular/human_txome_AdV41-A549-48h-1.GRCh38.polyA.tsv", "h48", ensembl)


# Summary checks (kept for parity)

for(df in list(h12, h24, h48, t12, t24, t48)) {
  print(summary(df$polya_length))
  print(getmode(round(df$polya_length, 0)))
}


# Combine all data

full_polya <- bind_rows(h12, h24, h48, t12, t24, t48) %>%
  mutate(
    contig = factor(contig),
    polya_length = as.numeric(polya_length)
  )


# Plot

p <- ggplot(full_polya, aes(
  x = factor(Run,
             levels = c("h12","h24","h48","t12","t24","t48"),
             labels = c(" 12 ", " 24 ", " 48 ", "12", "24", "48")),
  y = polya_length,
  fill = Run
)) +
  geom_violin(linewidth = 1) +
  geom_boxplot(width = 0.1, linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(values = c(
    h12 = "lightblue", h24 = "lightblue", h48 = "lightblue",
    t12 = "lightyellow", t24 = "lightyellow", t48 = "lightyellow"
  )) +
  ggtitle(plot_title) +
  xlab("hours post infection") +
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
ggsave(plot_file, plot = p, width = 15, height = 5, units = "in")
