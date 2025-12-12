library(biomaRt)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#setwd("")

# Connect once
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  version = 104
)

# Biotype filters to exclude
biotype_exclude <- c(
  "polymorphic_pseudogene", "artifact", "transcribed_unitary_pseudogene",
  "scaRNA", "snRNA", "Mt_tRNA", "ribozyme", "misc_RNA", "rRNA", "snoRNA",
  "non_stop_decay", "processed_transcript", "scRNA", "protein_coding_LoF",
  "processed_pseudogene", "unprocessed_pseudogene", "nonsense_mediated_decay",
  "retained_intron", "lncRNA"
)

library(tidyverse) ### Placed here to avoid conflicts with biomaRt package

process_sample <- function(df, sample_name, mart = ensembl) {
  
  df <- df %>%
    mutate(contig_no_version = sub("\\..*", "", contig))
  
  enst_ids <- unique(df$contig_no_version)
  
  annot <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                   "external_gene_name", "transcript_biotype"),
    filters = "ensembl_transcript_id",
    values = enst_ids,
    mart = mart
  )
  
  df <- df %>%
    left_join(annot, by = c("contig_no_version" = "ensembl_transcript_id")) %>%
    filter(qc_tag == "PASS") %>%
    filter(!grepl("Mt_\\|?$", transcript_biotype)) %>%
    filter(
      !is.na(transcript_biotype),
      !(transcript_biotype %in% biotype_exclude),
      !grepl("^MT-", external_gene_name, ignore.case = TRUE)
    ) %>%
    select(Transcript = contig_no_version, PolyA_Length = polya_length) %>%
    mutate(
      Transcript = sub("\\..*", "", Transcript),
      PolyA_Length = as.numeric(PolyA_Length),
      Run = sample_name
    )
  
  return(df)
}


# Load input data (ensure cellular and viral directories are in the same folder as this script)
A549  <- read.table("cellular/human_txome_A549-1.polyA.tsv", sep = "\t", header = TRUE)
MeWo  <- read.table("cellular/human_txome_MeWo-1.polyA.tsv", sep = "\t", header = TRUE)
NHDF  <- read.table("cellular/human_txome_NHDF-new.tsv", sep = "\t", header = TRUE)

# Process
A549_clean <- process_sample(A549, "A549")
MeWo_clean <- process_sample(MeWo, "MeWo")
NHDF_clean <- process_sample(NHDF, "NHDF")

# Combine
full_polya <- bind_rows(A549_clean, MeWo_clean, NHDF_clean)

# Plot
title <- "Uninfected human cell lines"
save  <- "Uninfected-human-polyA-protCodingonly.pdf"

p <- ggplot(full_polya, aes(
  x = factor(Run, levels = c("A549", "MeWo", "NHDF")),
  y = PolyA_Length,
  fill = Run
)) +
  geom_violin(lwd = 1) +
  geom_boxplot(width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(values = c("A549"="skyblue1", "MeWo"="skyblue1", "NHDF"="skyblue1")) +
  ggtitle(title) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    panel.grid.major = element_line(linewidth = 0.5, colour = "darkgrey"),
    panel.grid.minor = element_line(linewidth = 0.1, colour = "darkgrey"),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.text.x  = element_text(size = 25),
    axis.text.y  = element_text(size = 20),
    axis.title   = element_text(size = 25)
  )

print(p)
ggsave(save, plot = p, width = 8, height = 5, units = "in")
