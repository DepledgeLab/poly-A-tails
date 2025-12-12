library(tidyverse)

#setwd("")

#functions
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

clean_polyA_df <- function(df) {
  df %>%
    filter(qc_tag == "PASS") %>%
    select(contig, polya_length) %>%
    rename(Transcript = contig, PolyA_Length = polya_length) %>%
    mutate(
      Transcript = gsub("\\..*", "", Transcript),
      PolyA_Length = as.numeric(PolyA_Length)
    )
}

#Load and process PAN

pan_raw <- read.table("fig3_data/matched.pan.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(pan_raw) <- c(
  "readname","contig","position","leader_start","adapter_start",
  "polya_start","transcript_start","read_rate","polya_length","qc_tag"
)

pan <- pan_raw %>%
  clean_polyA_df() %>%
  mutate(Run = "PAN")


#Load non-PAN and split Kaposin
np_raw <- read.table("fig3_data/matched.not_pan.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(np_raw) <- c(
  "readname","contig","position","leader_start","adapter_start",
  "polya_start","transcript_start","read_rate","polya_length","qc_tag"
)

# Kaposin region definition
kap_mask <- np_raw$position >= 117500 & np_raw$position <= 117600

kaposin <- np_raw[kap_mask, ] %>% clean_polyA_df() %>% mutate(Run = "Kaposin")

non_pan_no_kap <- np_raw[!kap_mask, ] %>%
  clean_polyA_df() %>%
  mutate(Run = "viral mRNAs")


### 4. Load MALAT and NEAT ncRNAs


nc_raw <- read.table("fig3_data/KSHV-iSLK-72h-1-hac.v6.5.7.MALAT_NEAT.polyA.tsv",
                     sep="\t", header=FALSE, stringsAsFactors=FALSE)

colnames(nc_raw) <- c(
  "readname","contig","position","leader_start","adapter_start",
  "polya_start","transcript_start","read_rate","polya_length","qc_tag"
)

# Keep contig until after filtering
nc_raw_pass <- nc_raw %>% filter(qc_tag == "PASS")

MALAT <- nc_raw_pass %>%
  filter(grepl("MALAT", contig)) %>%
  clean_polyA_df() %>%
  mutate(Run = "MALAT1")

NEAT <- nc_raw_pass %>%
  filter(grepl("NEAT", contig)) %>%
  clean_polyA_df() %>%
  mutate(Run = "NEAT")


### 5. Combine all datasets


full_polya <- bind_rows(
  pan,
  kaposin,
  non_pan_no_kap,
  NEAT  # intentionally excluding MALAT to match your earlier choice
)

full_polya$Run <- factor(
  full_polya$Run,
  levels = c("PAN","Kaposin","viral mRNAs","NEAT"),
  labels = c("PAN","Kaposin","KSHV mRNAs","NEAT1")
)


### 6. Plot


title <- "KSHV, iSLK, 72hpi"

p <- ggplot(full_polya, aes(
  x = Run,
  y = PolyA_Length,
  fill = Run
)) +
  geom_violin(lwd = 1) +
  geom_boxplot(width = 0.1, lwd = 1) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  ylim(0, 1000) +
  scale_fill_manual(values = c(
    "PAN"="coral1",
    "Kaposin"="lightyellow",
    "KSHV mRNAs"="lightyellow",
    "NEAT1"="lightblue1"
  )) +
  ggtitle(title) +
  xlab(NULL) +
  ylab("poly(A) tail length") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(linewidth = 0.5, color = "darkgrey"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "darkgrey"),
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(color = "black", size = 1),
    axis.line.y = element_line(color = "black", size = 1),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 20)
  )

print(p)

ggsave("KSHV-PAN-Kaposin-analysis.pdf", p, width = 10, height = 5, units = "in")
