library(tidyverse)

#setwd("")

### Parameters
transcript <- "RL2-1"
title_gene <- "ICP47 (US12)"
save_gene <- "HSV1-ICP47-timecourse-polyA.pdf"

title_all <- "HSV-1 mRNAs"
save_all <- "HSV1-all-mRNAs-timecourse-polyA.pdf"


### Generic cleaning function for polyA files
clean_polyA <- function(df) {
  df %>%
    filter(qc_tag == "PASS") %>%
    select(contig, polya_length) %>%
    rename(Transcript = contig, PolyA_Length = polya_length) %>%
    mutate(
      Transcript = gsub("\\..*", "", Transcript),
      PolyA_Length = as.numeric(PolyA_Length)
    )
}



### Load the three datasets (3h, 6h, 12h)
files <- list(
  `3`  = "HSVtx/KOS_3-1.KOStx.polyA.tsv",
  `6`  = "HSVtx/Kos_6-1.KOStx.polyA.tsv",
  `12` = "HSVtx/Kos_12-1.KOStx.polyA.tsv"
)

raw_list <- map(files, ~ read.table(.x, sep="\t", header=TRUE, stringsAsFactors=FALSE))

### Process transcript of interest (ICP47 / RL2-1)
gene_list <- map2(
  raw_list,
  names(raw_list),
  ~ {
    df <- .x %>% filter(qc_tag == "PASS")
    df_gene <- df %>% filter(grepl(transcript, contig))
    clean_polyA(df_gene) %>% mutate(Time = .y)
  }
)

gene_df <- bind_rows(gene_list)

### Process all viral reads (for second plot)
all_list <- map2(
  raw_list,
  names(raw_list),
  ~ clean_polyA(.x) %>% mutate(Time = .y)
)

all_df <- bind_rows(all_list)

### Plot 1: transcript of interest (ICP47)
p <- ggplot(gene_df,
            aes(x = factor(Time, levels=c("3","6","12")),
                y = PolyA_Length,
                fill = Time)) +
  geom_violin(lwd=1) +
  geom_boxplot(width=0.1, lwd=1) +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  scale_fill_manual(values=c("3"="lightyellow","6"="lightyellow","12"="lightyellow")) +
  ylim(0, 1000) +
  ggtitle(title_gene) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 25),
    axis.text.x  = element_text(size = 25),
    axis.text.y  = element_text(size = 20),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "darkgrey"),
    panel.grid.minor = element_line(color = "darkgrey", linewidth = 0.1),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black")
  )

ggsave(save_gene, plot=p, width=7.5, height=5, units="in")
print(p)



### Plot 2: all viral transcripts
q <- ggplot(all_df,
            aes(x = factor(Time, levels=c("3","6","12")),
                y = PolyA_Length,
                fill = Time)) +
  geom_violin(lwd=1) +
  geom_boxplot(width=0.1, lwd=1) +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  scale_fill_manual(values=c("3"="lightyellow","6"="lightyellow","12"="lightyellow")) +
  ylim(0, 1000) +
  ggtitle(title_all) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.y = element_text(size = 25),
    axis.text.x  = element_text(size = 25),
    axis.text.y  = element_text(size = 20),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "darkgrey"),
    panel.grid.minor = element_line(color = "darkgrey", linewidth = 0.1),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(colour = "black")
  )

ggsave(save_all, plot=q, width=7.5, height=5, units="in")
print(q)
