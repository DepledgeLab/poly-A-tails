library(tidyverse)
library(biomaRt)



# Mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Paths
#setwd("")

# Load input data (ensure cellular and viral directories are in the same folder as this script)
human_df <- read.table(
  "cellular/human_txome_EBV.polyA.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

virus_df <- read.table(
  "viral/viral_EBV.polyA.tsv",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Edit as needed 
plot_title <- "EBV-B-95-8-BAC, HEK293, 48 hpr"
plot_save  <- "EBV-human-polyA-protCodingonly.pdf"

# Ensembl
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# Remove transcript version from human table
human_df <- human_df %>%
  mutate(contig_no_version = sub("\\..*", "", contig))

# Pull annotation
annot <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "external_gene_name",
    "transcript_biotype"
  ),
  filters = "ensembl_transcript_id",
  values = unique(human_df$contig_no_version),
  mart = ensembl
)

# Join annotation
merged <- human_df %>%
  dplyr::left_join(
    annot,
    by = c("contig_no_version" = "ensembl_transcript_id")
  )

# Extract and filter
filtered <- merged %>%
  select(
    contig_no_version,
    polya_length,
    qc_tag,
    external_gene_name,
    transcript_biotype
  ) %>%
  filter(qc_tag == "PASS") %>%
  filter(!grepl("Mt_\\|?$", transcript_biotype)) %>%
  filter(
    !is.na(transcript_biotype),
    !grepl(
      paste(
        c(
          "polymorphic_pseudogene", "artifact", "transcribed_unitary_pseudogene",
          "scaRNA", "snRNA", "Mt_tRNA", "ribozyme", "misc_RNA", "rRNA", "snoRNA",
          "non_stop_decay", "processed_transcript", "scRNA", "protein_coding_LoF",
          "processed_pseudogene", "unprocessed_pseudogene",
          "nonsense_mediated_decay", "retained_intron", "lncRNA"
        ),
        collapse = "|"
      ),
      transcript_biotype
    ),
    !grepl("^MT-", external_gene_name, ignore.case = TRUE)
  )

# Final clean formatting
human_clean <- filtered %>%
  select(Transcript = contig_no_version, PolyA_Length = polya_length) %>%
  mutate(
    Transcript = sub("\\..*", "", Transcript),
    PolyA_Length = as.numeric(PolyA_Length),
    Run = "Human"
  )

# Virus processing
virus_clean <- virus_df %>%
  filter(qc_tag == "PASS") %>%
  select(-c(1, 3:8, 10)) %>%
  rename(Transcript = contig, PolyA_Length = polya_length) %>%
  mutate(
    Transcript = sub("\\..*", "", Transcript),
    PolyA_Length = as.numeric(PolyA_Length),
    Run = "Viral"
  )

# Merge both for plotting
full_polya <- bind_rows(human_clean, virus_clean)

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
    values = c(Human = "lightblue", Viral = "lightyellow")
  ) +
  ggtitle(plot_title) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20, colour = "black"),
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
ggsave(plot_save, plot = p, width = 5, height = 5, units = "in")
