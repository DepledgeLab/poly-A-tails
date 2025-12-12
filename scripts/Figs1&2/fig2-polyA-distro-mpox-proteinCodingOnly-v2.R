library(tidyverse)
library(biomaRt)
library(stringr)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#setwd("")

plot_title <- "Mpox, NHDF, MOI=1"
plot_file  <- "Mpox-timecourse-polyA-protCodingonly.pdf"

ensembl <- useEnsembl(
  biomart  = "ensembl",
  dataset  = "hsapiens_gene_ensembl",
  version  = 104
)

#process viral files
process_virus <- function(path, label) {
  read_tsv(path, show_col_types = FALSE) %>%
    dplyr::filter(qc_tag == "PASS") %>%
    dplyr::select(contig, polya_length) %>%
    dplyr::mutate(
      contig = sub("\\..*", "", contig),
      polya_length = as.numeric(polya_length),
      Run = label
    )
}

#process human files
process_human <- function(path, label, mart) {
  df <- read_tsv(path, show_col_types = FALSE) %>%
    dplyr::mutate(contig_no_version = sub("\\..*", "", contig))
  
  enst_ids <- unique(df$contig_no_version)
  
  annot <- getBM(
    attributes = c(
      "ensembl_transcript_id",
      "ensembl_gene_id",
      "external_gene_name",
      "transcript_biotype"
    ),
    filters = "ensembl_transcript_id",
    values = enst_ids,
    mart = mart
  ) %>%
    as_tibble()
  
  # biotypes to exclude (kept as a regex for one-line filtering)
  exclude_re <- paste0(
    "polymorphic_pseudogene|artifact|transcribed_unitary_pseudogene|",
    "scaRNA|snRNA|Mt_tRNA|ribozyme|misc_RNA|rRNA|snoRNA|",
    "non_stop_decay|processed_transcript|scRNA|protein_coding_LoF|",
    "processed_pseudogene|unprocessed_pseudogene|nonsense_mediated_decay|",
    "retained_intron|lncRNA"
  )
  
  df %>%
    dplyr::left_join(annot, by = c("contig_no_version" = "ensembl_transcript_id")) %>%
    dplyr::filter(qc_tag == "PASS") %>%
    dplyr::filter(!is.na(transcript_biotype)) %>%
    dplyr::filter(!str_detect(transcript_biotype, exclude_re)) %>%
    dplyr::filter(!str_detect(external_gene_name, regex("^MT-", ignore_case = TRUE))) %>%
    dplyr::select(contig = contig_no_version, polya_length) %>%
    dplyr::mutate(
      contig = sub("\\..*", "", contig),
      polya_length = as.numeric(polya_length),
      Run = label
    )
}

# # Load input data (ensure viral directory is in the same folder as this script)
t12 <- process_virus("viral/viral_MPX-NL1-NHDF-CTRL-4h.MPXV-UK.polyA.tsv", "t12")
t24 <- process_virus("viral/viral_MPX-NL1-NHDF-CTRL-10h.MPXV-UK.polyA.tsv", "t24")

# # Load input data (ensure cellular directory is in the same folder as this script)
h12 <- process_human("cellular/human_txome_MPX-NL1-NHDF-CTRL-4h.GRCh38.polyA.tsv", "h12", ensembl)
h24 <- process_human("cellular/human_txome_MPX-NL1-NHDF-CTRL-10h.GRCh38.polyA.tsv", "h24", ensembl)

#quick summaries (parity with original)
list(h12 = h12, h24 = h24, t12 = t12, t24 = t24) %>%
  imap(~{
    message("Summary for: ", .y)
    print(summary(.x$polya_length))
    message("Mode (rounded): ", getmode(round(.x$polya_length, 0)))
    invisible(NULL)
  })

#combine for plotting
full_polya <- bind_rows(h12, h24, t12, t24) %>%
  mutate(
    contig = factor(contig),
    PolyA_Length = polya_length  # keep naming consistent with your plotting code
  )

#plot
p <- ggplot(full_polya, aes(
  x = factor(Run, levels = c("h12","h24","t12","t24"),
             labels = c(" 4 ", " 10 ", "4", "10")),
  y = PolyA_Length,
  fill = Run
)) +
  geom_violin(linewidth = 1) +
  geom_boxplot(width = 0.1, linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(values = c(
    "h12" = "lightblue", "h24" = "lightblue",
    "t12" = "lightyellow", "t24" = "lightyellow"
  )) +
  ggtitle(plot_title) +
  xlab("hours post infection") +
  ylab("poly(A) tail length") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )

print(p)
ggsave(plot_file, plot = p, width = 10, height = 5, units = "in")
