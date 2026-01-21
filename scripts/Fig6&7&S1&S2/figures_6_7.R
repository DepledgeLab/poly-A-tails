###############################################################################\
#                           Figure 6/7 generation (& S1/2)                  ####
###############################################################################/



###############################################################################\
## Notes on running this script ####
#
# For performance reasons, this scripts was run on the authors universities HPC,
# via rocker/ RStudio Server. More detailed information regarding the used
# container and conda environments can be found in the associated files in the
# folder containing this script.
# General installation information for Ninetails can be found on their GitHub
# page (https://github.com/LRB-IIMCB/ninetails).

## Further notes:
# - In the context of this script, 'decorated' poly(A) tails are synonymously
# used with 'mixed' poly(A) tails.
# - Some individual figures from this script have been joined to a single figure
# in the final paper.
# - Column names of the GLM output have been renamed in the supplementary files
# for comprehensibility.


###############################################################################/



###############################################################################\
## Libraries & Env                                                          ####
###############################################################################/

rm(list = ls())

library(tidyr)
library(ggplot2)
library(viridis)
library(rstatix)
library(ggpubr)
library(stringr)
library(forcats)
library(ggthemes)
library(FSA)
library(ggsignif)
library(cowplot)
library(patchwork)
library(grDevices)
library(zoo)
library(dplyr)
library(glmmTMB)

# Set env based on your specific paths and environment names
setwd("/path/to/poly-A-tails/scripts/Fig6&7&S1&S2")
Sys.setenv(RETICULATE_CONDA="/path/to/Anaconda3/2022.05/bin/conda")
reticulate::use_condaenv("r-ninetails")


###############################################################################/



###############################################################################\
## Functions                                                                ####
###############################################################################/

# Fit NB GLMM per transcript
fit_one_transcript <- function(tid) {
  dat_t <- dplyr::bind_rows(
    baseline_df %>% mutate(ensembl_transcript_id_full = tid),
    latest_df %>% filter(ensembl_transcript_id_full == tid) %>%
      select(ensembl_transcript_id_full, sample, y, E, n_reads, comp)
  )
  
  # Remove any rows with missing values
  dat_t <- dat_t %>% tidyr::drop_na(y, E, sample, comp)
  dat_t <<- dat_t
  
  # Sanity checks
  if (!all(c(0,1) %in% dat_t$comp))
    warning(sprintf("Transcript %s: comp missing a level (levels: %s)",
                    tid,
                    paste(unique(dat_t$comp),collapse=",")))
  if (!("sample" %in% names(dat_t)))
    stop("Column 'sample' not found in dat_t")
  if (sum(dat_t$y) == 0)
    warning(sprintf("Transcript %s: all counts are zero; NB GLMM may fail",
                    tid))
  if (any(!is.finite(dat_t$E)) || any(dat_t$E <= 0))
    stop(sprintf("Transcript %s: non-finite or non-positive exposure E", tid))
  
  dat_t$sample <- factor(dat_t$sample)
  
  # Fit NB GLMM with offset for exposure
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(y ~ comp + (1 | sample),
              offset = log(E),
              family = nbinom2,
              data = dat_t),
      warning = function(w) {
        message(sprintf("Warning for transcript %s: %s",
                        tid,
                        conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    i <<- i + 1
    warning(sprintf("NB GLMM fit failed for transcript %s; returning NA results and flagging as 'glmmTMB_nbinom2_failed'.", tid))
    return(
      data.frame(
        transcript_id = tid,
        method = "glmmTMB_nbinom2_failed",
        rate_ratio = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_,
        total_reads_latest = sum(dat_t$n_reads[dat_t$comp == 1]),
        total_E_latest = sum(dat_t$E[dat_t$comp == 1]),
        total_y_latest = sum(dat_t$y[dat_t$comp == 1]),
        total_reads_mock = sum(dat_t$n_reads[dat_t$comp == 0]),
        total_E_mock = sum(dat_t$E[dat_t$comp == 0]),
        total_y_mock = sum(dat_t$y[dat_t$comp == 0])
      )
    )
  } else {
    summ <- summary(fit)
    beta <- summ$coefficients$cond["comp", "Estimate"]
    pval <- summ$coefficients$cond["comp", "Pr(>|z|)"]
    ci   <- suppressWarnings(confint(fit, parm = "comp", level = 0.95))
    
    # handle possible matrix/vector return from confint
    if (is.matrix(ci)) {
      lcl <- ci[1, 1]; ucl <- ci[1, 2]
    } else {
      lcl <- ci[1]; ucl <- ci[2]
    }
    return(
      data.frame(
        transcript_id = tid,
        method = "glmmTMB_nbinom2",
        rate_ratio = exp(beta),
        conf_low = exp(lcl),
        conf_high = exp(ucl),
        p_value = pval,
        total_reads_latest = sum(dat_t$n_reads[dat_t$comp == 1]),
        total_E_latest = sum(dat_t$E[dat_t$comp == 1]),
        total_y_latest = sum(dat_t$y[dat_t$comp == 1]),
        total_reads_mock = sum(dat_t$n_reads[dat_t$comp == 0]),
        total_E_mock = sum(dat_t$E[dat_t$comp == 0]),
        total_y_mock = sum(dat_t$y[dat_t$comp == 0])
      )
    )
  }
}


# Get variable name as string
get_var_name <- function(var) {
  deparse(substitute(var))
}


# Apply GLM stats with data wrangling on a specific nucleotide type
apply_glm <- function(merged_data_reclass_hcmv_human,
                      to_test_timepopint,
                      nucleotide_pred){
  
  nucleotide_pred_col <- get_var_name(nucleotide_pred)
  
  # Check required columns
  stopifnot(all(c("sample",
                  "timepoint",
                  "group",
                  "ensembl_transcript_id_full",
                  "polya_round",
                  nucleotide_pred_col) %in%
                  names(merged_data_reclass_hcmv_human)))
  
  # Aggregate pooled human Mock baseline (pooled across all human transcripts)
  baseline_df <- merged_data_reclass_hcmv_human %>%
    filter(timepoint == "Mock", group == "human") %>%
    group_by(sample, group) %>%
    summarise(y = sum(nucleotide_pred),
              E = sum(polya_round),
              n_reads = n(),
              .groups = "drop") %>%
    mutate(comp = 0L)
  
  if (nrow(baseline_df) == 0)
    stop("No human Mock baseline rows found.")
  
  # Aggregate latest-timepoint counts per transcript x sample
  latest_df <- merged_data_reclass_hcmv_human %>%
    filter(timepoint == to_test_timepopint) %>%
    group_by(ensembl_transcript_id_full, sample, group) %>%
    summarise(y = sum(nucleotide_pred),
              E = sum(polya_round),
              n_reads = n(),
              .groups = "drop") %>%
    mutate(comp = 1L)
  
  if (nrow(latest_df) == 0)
    stop("No rows found for latest timepoint: ", to_test_timepopint)
  
  # Apply depth filter:
  # require minimum total latest reads per transcript across replicates
  transcripts_to_test <- latest_df %>%
    group_by(ensembl_transcript_id_full) %>%
    summarise(total_reads_latest = sum(n_reads),
              total_E_latest = sum(E),
              .groups = "drop") %>%
    filter(total_reads_latest >= min_number_reads) %>%
    pull(ensembl_transcript_id_full)
  
  if (length(transcripts_to_test) == 0)
    stop("No transcripts pass the depth threshold of ",
         min_number_reads,
         " total latest reads.")
  
  latest_df <- latest_df %>%
    filter(ensembl_transcript_id_full %in% transcripts_to_test)
  
  # Ensure exposure > 0
  baseline_df <- baseline_df %>%
    filter(E > 0)
  latest_df   <- latest_df   %>%
    filter(E > 0)
  
  # Benjamini-Hochberg correction
  results_list <- lapply(transcripts_to_test, fit_one_transcript)
  results_df <- bind_rows(results_list) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    arrange(q_value, p_value)
  
  return(list(results_list, results_df))
  
  
}

# Annotate GML reulsts dataframe
annot_glm_df <- function(results_df, nucleotide) {
  
  results_df_2 <- results_df %>%
    left_join(group_map, by = "transcript_id") %>% 
    mutate(transcript_name = case_when(
      group == "human" ~ transcript_id,
      group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", transcript_id),
    )) %>% 
    mutate(dec_type = nucleotide)
  
  return(results_df_2)
  
}


# Get significant GLM results from dataframe
get_GLM_sign <- function(results_df, nucleotide) {
  
  results_df_sig <- results_df %>%
    filter(q_value <= 0.05) %>% # alpha threshold
    left_join(group_map, by = "transcript_id") %>% 
    mutate(transcript_name = case_when(
      group == "human" ~ transcript_id,
      group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", transcript_id),
    )) %>% 
    mutate(dec_type = nucleotide)
  
  return(results_df_sig)
  
}

# Change facet abels for the Manhatten-like plots
change_facet_labels <- function(variable, value){
  return(label_list[value])
}

###############################################################################/



###############################################################################\
## Base Paths and general Vars                                              ####
###############################################################################/

### Set Path ###################################################################

POLYA_path <- "./data/"


### Set plot themes and aesthetics #############################################

# Plot theme
my_theme <- theme(
  axis.title.x = element_text(size = 25),
  axis.text.x  = element_text(size = 25, color = "black"),
  axis.line.x  = element_line(color = "black"),
  axis.ticks.x = element_line(color = "black"),
  
  axis.title.y = element_text(size = 25),
  axis.text.y  = element_text(size = 25, color = "black"),
  axis.line.y  = element_line(color = "black"),
  axis.ticks.y = element_line(color = "black"),
  
  axis.ticks.length = unit(5, "pt"),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
  panel.grid.major.x = element_blank(),
  
  legend.text = element_text(size = 15, color = "black"),
  legend.title = element_text(size = 20, color = "black"),
  legend.background = element_blank(),
  
  strip.text = element_text(size = 20),
  strip.text.x = element_text(size = 20),
  strip.text.y = element_text(size = 20),
  strip.background = element_blank(),
  text = element_text(size = 20),
  strip.placement = "outside"
)

#textbasesize <- 1 # fits for ~ 3.8/2.4 inch
#textbasesize <- 11 # fits for 1600:900
textbasesize <- 13 # fits for 2400:4800

my_theme_rel <- theme(
  text = element_text(size = textbasesize, color = "black"),
  axis.title = element_text(size = textbasesize, color = "black"),
  axis.text = element_text(size = textbasesize, color = "black"),
  strip.text = element_text(size = textbasesize, color = "black"),
  axis.line = element_line(linewidth = textbasesize, color = "black"),
  axis.ticks = element_line(linewidth = textbasesize, color = "black"),
  
  axis.title.x = element_text(size = rel(2.5)),
  axis.text.x  = element_text(size = rel(2)),
  axis.line.x  = element_line(linewidth = rel(0.125)),
  axis.ticks.x = element_line(linewidth = rel(0.125)),
  
  axis.title.y = element_text(size = rel(2.5)),
  axis.text.y  = element_text(size = rel(2)),
  axis.line.y  = element_line(linewidth = rel(0.125)),
  axis.ticks.y = element_line(linewidth = rel(0.125)),
  
  plot.title = element_text(size = rel(2.75), hjust = 0.5),
  
  axis.ticks.length = unit(1 * textbasesize, "pt"),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.grid.major.y = element_line(color = "grey90", linewidth = rel(0.125)),
  panel.grid.major.x = element_blank(),
  
  legend.text = element_text(size = rel(2)),
  legend.title = element_text(size = rel(2), face = "bold"),
  legend.background = element_blank(),
  
  strip.text.x = element_text(size = rel(2.25)),
  strip.text.y = element_text(size = rel(2.25)),
  strip.text.y.right = element_text(size = rel(1)),
  strip.background = element_blank(),
  strip.placement = "outside"
)

my_theme_rel_pie <- theme(
  text = element_text(size = textbasesize, color = "black"),
  strip.text = element_text(size = textbasesize, color = "black"),
  
  plot.title = element_text(size = rel(2.75), hjust = 0.5),
  
  legend.background = element_blank(),
  legend.title = element_text(size = rel(2), face = "bold"),
  legend.text = element_text(size = rel(2)),
  
  strip.text.x = element_text(size = rel(2.25)),
  strip.background = element_blank()
)

# label sizes
labelsize <- 2 * textbasesize
boxplot_pointsize <- 0.65 * textbasesize
boxplot_linesize <- 0.25 * textbasesize
boxplot_sigsize <- 1.65 * textbasesize
manh_pointsize <- 1.2 * textbasesize
manh_sigsize <- 2.1 * textbasesize
manh_guidesize <- 2 * textbasesize
manh_linesize <- 0.2 * textbasesize
plotlabelsize <- 3 * textbasesize
lineplot_linesize <- 0.1 * textbasesize


## Colors
# Viral
color_vir_24 <- "#FF00F5"
color_vir_48 <- "#B300AC"
color_vir_72 <- "#3E0056"
color_vir_72_si <- "#50006F"

# Human
color_hum_M <- "#00DDFF"
color_hum_24 <- "#0081D6"
color_hum_48 <- "#0024B9"
color_hum_72 <- "#0C005F"
color_hum_72_si <- "#070079"

# Others
color_2_7 <- "coral1"
color_UL69 <- "goldenrod3"
color_UL22A <- "goldenrod1"
color_UL132 <- "goldenrod3"
color_UL77 <- "goldenrod4"
color_B2M <- "#87C7FF"
color_HSP <- "#0024B9"

color_dec <- "#47098E"
color_dec2 <- "#24004A"
color_bla_hum <- "#87C7FF"
color_bla_hum2 <- "#2268A7"
color_bla_vir <- "#FFEF71"
color_bla_vir2 <- "#AD9E2C"

###############################################################################\
## Loading in data                                                          ####
###############################################################################/


# TB40 POLYA (HCMV)
TB40_24h_class <- read.table(paste0(POLYA_path,
                                    "TB40_24h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.class_data_reclass.tsv"),
                             sep = "\t",
                             header = TRUE)
TB40_48h_class <- read.table(paste0(POLYA_path,
                                    "TB40_48h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.class_data_reclass.tsv"),
                             sep = "\t",
                             header = TRUE)
TB40_72h_class <- read.table(paste0(POLYA_path,
                                    "TB40_72h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.class_data_reclass.tsv"),
                             sep = "\t",
                             header = TRUE)

TB40_24h_residue <- read.table(paste0(POLYA_path,
                                      "TB40_24h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.residue_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)
TB40_48h_residue <- read.table(paste0(POLYA_path,
                                      "TB40_48h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.residue_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)
TB40_72h_residue <- read.table(paste0(POLYA_path,
                                      "TB40_72h.TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF.residue_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)

# TB40 siCTRL (HCMV)
TB40CTRL_72h_1_class <- read.table(paste0(POLYA_path,
                                          "TB40_72h_CTRL-1.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf.human-noMT_protEnc.hcmv.ninetails_pF.class_data_reclass.tsv"),
                                   sep = "\t",
                                   header = TRUE)
TB40CTRL_72h_2_class <- read.table(paste0(POLYA_path,
                                          "TB40_72h_CTRL-2.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf.human-noMT_protEnc.hcmv.ninetails_pF.class_data_reclass.tsv"),
                                   sep = "\t",
                                   header = TRUE)

TB40CTRL_72h_1_residue <- read.table(paste0(POLYA_path,
                                            "TB40_72h_CTRL-1.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf.human-noMT_protEnc.hcmv.ninetails_pF.residue_data_reclass.tsv"),
                                     sep = "\t",
                                     header = TRUE)
TB40CTRL_72h_2_residue <- read.table(paste0(POLYA_path,
                                            "TB40_72h_CTRL-2.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf.human-noMT_protEnc.hcmv.ninetails_pF.residue_data_reclass.tsv"),
                                     sep = "\t",
                                     header = TRUE)



# KSHV-iSLK-72h-1
KSHV_72h_1_class <- read.table(paste0(POLYA_path,
                                      "KSHV-iSLK-72h-1.guppy.hac.6.1.7.KSHVtxome.uf.ninetails_pF.class_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)

KSHV_72h_1_residue <- read.table(paste0(POLYA_path,
                                        "KSHV-iSLK-72h-1.guppy.hac.6.1.7.KSHVtxome.uf.ninetails_pF.residue_data_reclass.tsv"),
                                 sep = "\t",
                                 header = TRUE)


# NHDFpolyA-RN7SK-002 (Mock)
NHDFRN7SK_class <- read.table(paste0(POLYA_path,
                                     "NHDFpolyA-RN7SK-002.guppy.hac.6.1.7gencode.v47_ENO2.tx.uf.polyA.human-noMT_protEnc.ENO2.RN7SK.ninetails_pF.class_data_reclass.tsv"),
                              sep = "\t",
                              header = TRUE)

NHDFRN7SK_residue <- read.table(paste0(POLYA_path,
                                       "NHDFpolyA-RN7SK-002.guppy.hac.6.1.7gencode.v47_ENO2.tx.uf.polyA.human-noMT_protEnc.ENO2.RN7SK.ninetails_pF.residue_data_reclass.tsv"),
                                sep = "\t",
                                header = TRUE)


# IVT RN7SKpolyA-002
IVTRN7SK_class <- read.table(paste0(IVTRN7SK_path,
                                    "RN7SKpolyA-002.guppy.hac.6.1.7.RN7SK.tx.uf.ninetails_pF.class_data_reclass.tsv"),
                             sep = "\t",
                             header = TRUE)

IVTRN7SK_residue <- read.table(paste0(IVTRN7SK_path,
                                      "RN7SKpolyA-002.guppy.hac.6.1.7.RN7SK.tx.uf.ninetails_pF.residue_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)

# HSV2-ARPE19-10h-1
HSV2_10h_class <- read.table(paste0(POLYA_path,
                                    "HSV2-ARPE19-10h-1.guppy.hac.6.1.7.HSV2-MS.tx.uf.ninetails_pF.class_data_reclass.tsv"),
                             sep = "\t",
                             header = TRUE)

HSV2_10h_residue <- read.table(paste0(POLYA_path,
                                      "HSV2-ARPE19-10h-1.guppy.hac.6.1.7.HSV2-MS.tx.uf.ninetails_pF.residue_data_reclass.tsv"),
                               sep = "\t",
                               header = TRUE)

# EMC1-MeWo-96h-polyA (VZV)
Dumas_96h_class <- read.table(paste0(POLYA_path,
                                     "EMC1-MeWo-96h-polyA.guppy.hac.6.1.7.dumas.tx.uf.ninetails_pF.class_data_reclass.tsv"),
                              sep = "\t",
                              header = TRUE)
Dumas_96h_residue <- read.table(paste0(POLYA_path,
                                       "EMC1-MeWo-96h-polyA.guppy.hac.6.1.7.dumas.tx.uf.ninetails_pF.residue_data_reclass.tsv"),
                                sep = "\t",
                                header = TRUE)



###############################################################################\
## Merging the data to single dataframes                                    ####
###############################################################################/

# Merge class and residue dataframes

merged_TB40_reclass_24h <- ninetails::merge_nonA_tables(class_data=TB40_24h_class,
                                                        residue_data=TB40_24h_residue,
                                                        pass_only=FALSE)
merged_TB40_reclass_48h <- ninetails::merge_nonA_tables(class_data=TB40_48h_class,
                                                        residue_data=TB40_48h_residue,
                                                        pass_only=FALSE)
merged_TB40_reclass_72h <- ninetails::merge_nonA_tables(class_data=TB40_72h_class,
                                                        residue_data=TB40_72h_residue,
                                                        pass_only=FALSE)
merged_KSHV_reclass_72h_1 <- ninetails::merge_nonA_tables(class_data=KSHV_72h_1_class,
                                                          residue_data=KSHV_72h_1_residue,
                                                          pass_only=FALSE)
merged_NHDFRN7SK_reclass <- ninetails::merge_nonA_tables(class_data=NHDFRN7SK_class,
                                                         residue_data=NHDFRN7SK_residue,
                                                         pass_only=FALSE)
merged_IVTRN7SK_reclass <- ninetails::merge_nonA_tables(class_data=IVTRN7SK_class,
                                                        residue_data=IVTRN7SK_residue,
                                                        pass_only=FALSE)
merged_TB40CTRL_reclass_72h_1 <- ninetails::merge_nonA_tables(class_data=TB40CTRL_72h_1_class,
                                                              residue_data=TB40CTRL_72h_1_residue,
                                                              pass_only=FALSE)
merged_TB40CTRL_reclass_72h_2 <- ninetails::merge_nonA_tables(class_data=TB40CTRL_72h_2_class,
                                                              residue_data=TB40CTRL_72h_2_residue,
                                                              pass_only=FALSE)
merged_HSV2_reclass_10h <- ninetails::merge_nonA_tables(class_data=HSV2_10h_class,
                                                        residue_data=HSV2_10h_residue,
                                                        pass_only=FALSE)
merged_Dumas_reclass_96h <- ninetails::merge_nonA_tables(class_data=Dumas_96h_class,
                                                         residue_data=Dumas_96h_residue,
                                                         pass_only=FALSE)


# Summarise dataframes (per transcript)

merged_TB40_reclass_summarized_24h <- ninetails::summarize_nonA(merged_nonA_tables=merged_TB40_reclass_24h,
                                                                summary_factors="group",
                                                                transcript_id_column="ensembl_transcript_id_full")
merged_TB40_reclass_summarized_48h <- ninetails::summarize_nonA(merged_nonA_tables=merged_TB40_reclass_48h,
                                                                summary_factors="group",
                                                                transcript_id_column="ensembl_transcript_id_full")
merged_TB40_reclass_summarized_72h <- ninetails::summarize_nonA(merged_nonA_tables=merged_TB40_reclass_72h,
                                                                summary_factors="group",
                                                                transcript_id_column="ensembl_transcript_id_full")
merged_KSHV_reclass_summarized_72h_1 <- ninetails::summarize_nonA(merged_nonA_tables=merged_KSHV_reclass_72h_1,
                                                                  summary_factors="group",
                                                                  transcript_id_column="ensembl_transcript_id_full")
merged_NHDFRN7SK_reclass_summarized <- ninetails::summarize_nonA(merged_nonA_tables=merged_NHDFRN7SK_reclass,
                                                                 summary_factors="group",
                                                                 transcript_id_column="ensembl_transcript_id_full")
merged_IVTRN7SK_reclass_summarized <- ninetails::summarize_nonA(merged_nonA_tables=merged_IVTRN7SK_reclass,
                                                                summary_factors="group",
                                                                transcript_id_column="ensembl_transcript_id_full")
merged_TB40CTRL_reclass_summarized_72h_1 <- ninetails::summarize_nonA(merged_nonA_tables=merged_TB40CTRL_reclass_72h_1,
                                                                      summary_factors="group",
                                                                      transcript_id_column="ensembl_transcript_id_full")
merged_TB40CTRL_reclass_summarized_72h_2 <- ninetails::summarize_nonA(merged_nonA_tables=merged_TB40CTRL_reclass_72h_2,
                                                                      summary_factors="group",
                                                                      transcript_id_column="ensembl_transcript_id_full")
merged_HSV2_reclass_summarized_10h <- ninetails::summarize_nonA(merged_nonA_tables=merged_HSV2_reclass_10h,
                                                                summary_factors="group",
                                                                transcript_id_column="ensembl_transcript_id_full")
merged_Dumas_reclass_summarized_96h <- ninetails::summarize_nonA(merged_nonA_tables=merged_Dumas_reclass_96h,
                                                                 summary_factors="group",
                                                                 transcript_id_column="ensembl_transcript_id_full")

# Get decoration rates and frequencies

merged_TB40_reclass_summarized_24h$decoration_rate <-
  merged_TB40_reclass_summarized_24h$counts_nonA / merged_TB40_reclass_summarized_24h$counts_total
merged_TB40_reclass_summarized_48h$decoration_rate <-
  merged_TB40_reclass_summarized_48h$counts_nonA / merged_TB40_reclass_summarized_48h$counts_total
merged_TB40_reclass_summarized_72h$decoration_rate <-
  merged_TB40_reclass_summarized_72h$counts_nonA / merged_TB40_reclass_summarized_72h$counts_total
merged_KSHV_reclass_summarized_72h_1$decoration_rate <-
  merged_KSHV_reclass_summarized_72h_1$counts_nonA / merged_KSHV_reclass_summarized_72h_1$counts_total
merged_NHDFRN7SK_reclass_summarized$decoration_rate <-
  merged_NHDFRN7SK_reclass_summarized$counts_nonA / merged_NHDFRN7SK_reclass_summarized$counts_total
merged_IVTRN7SK_reclass_summarized$decoration_rate <-
  merged_IVTRN7SK_reclass_summarized$counts_nonA / merged_IVTRN7SK_reclass_summarized$counts_total
merged_TB40CTRL_reclass_summarized_72h_1$decoration_rate <-
  merged_TB40CTRL_reclass_summarized_72h_1$counts_nonA / merged_TB40CTRL_reclass_summarized_72h_1$counts_total
merged_TB40CTRL_reclass_summarized_72h_2$decoration_rate <-
  merged_TB40CTRL_reclass_summarized_72h_2$counts_nonA / merged_TB40CTRL_reclass_summarized_72h_2$counts_total
merged_HSV2_reclass_summarized_10h$decoration_rate <-
  merged_HSV2_reclass_summarized_10h$counts_nonA / merged_HSV2_reclass_summarized_10h$counts_total
merged_Dumas_reclass_summarized_96h$decoration_rate <-
  merged_Dumas_reclass_summarized_96h$counts_nonA / merged_Dumas_reclass_summarized_96h$counts_total


# Add timepoints

merged_TB40_reclass_24h$timepoint <- "24h"
merged_TB40_reclass_48h$timepoint <- "48h"
merged_TB40_reclass_72h$timepoint <- "72h"
merged_KSHV_reclass_72h_1$timepoint <- "72h"
merged_NHDFRN7SK_reclass$timepoint <- "Mock"
merged_IVTRN7SK_reclass$timepoint <- "None"
merged_TB40CTRL_reclass_72h_1$timepoint <- "72h"
merged_TB40CTRL_reclass_72h_2$timepoint <- "72h"
merged_HSV2_reclass_10h$timepoint <-  "10h"
merged_Dumas_reclass_96h$timepoint <- "96h"


merged_TB40_reclass_summarized_24h$timepoint <- "24h"
merged_TB40_reclass_summarized_48h$timepoint <- "48h"
merged_TB40_reclass_summarized_72h$timepoint <- "72h"
merged_KSHV_reclass_summarized_72h_1$timepoint <- "72h"
merged_NHDFRN7SK_reclass_summarized$timepoint <- "Mock"
merged_IVTRN7SK_reclass_summarized$timepoint <- "None"
merged_TB40CTRL_reclass_summarized_72h_1$timepoint <- "72h"
merged_TB40CTRL_reclass_summarized_72h_2$timepoint <- "72h"
merged_HSV2_reclass_summarized_10h$timepoint <- "10h"
merged_Dumas_reclass_summarized_96h$timepoint <- "96h"


TB40_24h_residue$timepoint <- "24h"
TB40_48h_residue$timepoint <- "48h"
TB40_72h_residue$timepoint <- "72h"
TB40CTRL_72h_1_residue$timepoint <- "72h"
TB40CTRL_72h_2_residue$timepoint <- "72h"
KSHV_72h_1_residue$timepoint <- "72h"
NHDFRN7SK_residue$timepoint <- "Mock"
IVTRN7SK_residue$timepoint <- "None"
HSV2_10h_residue$timepoint <- "10h"
Dumas_96h_residue$timepoint <- "96h"


# Add sample numbers

merged_TB40_reclass_24h$sample <- "TB40_POLYA_24"
merged_TB40_reclass_48h$sample <- "TB40_POLYA_48"
merged_TB40_reclass_72h$sample <- "TB40_POLYA_72"
merged_KSHV_reclass_72h_1$sample <- "KSHV_iSLK"
merged_NHDFRN7SK_reclass$sample <- "NHDFRN7SK"
merged_IVTRN7SK_reclass$sample <- "IVTRN7SK"
merged_TB40CTRL_reclass_72h_1$sample <- "TB40_siCTRL_1"
merged_TB40CTRL_reclass_72h_2$sample <- "TB40_siCTRL_2"
merged_HSV2_reclass_10h$sample <-  "HSV2"
merged_Dumas_reclass_96h$sample <- "VZV"


merged_TB40_reclass_summarized_24h$sample <- "TB40_POLYA_24"
merged_TB40_reclass_summarized_48h$sample <- "TB40_POLYA_48"
merged_TB40_reclass_summarized_72h$sample <- "TB40_POLYA_72"
merged_KSHV_reclass_summarized_72h_1$sample <- "KSHV_iSLK"
merged_NHDFRN7SK_reclass_summarized$sample <- "NHDFRN7SK"
merged_IVTRN7SK_reclass_summarized$sample <- "IVTRN7SK"
merged_TB40CTRL_reclass_summarized_72h_1$sample <- "TB40_siCTRL_1"
merged_TB40CTRL_reclass_summarized_72h_2$sample <- "TB40_siCTRL_2"
merged_HSV2_reclass_summarized_10h$sample <-  "HSV2"
merged_Dumas_reclass_summarized_96h$sample <- "VZV"


TB40_24h_residue$sample <- "TB40_POLYA_24"
TB40_48h_residue$sample <- "TB40_POLYA_48"
TB40_72h_residue$sample <- "TB40_POLYA_72"
TB40CTRL_72h_1_residue$sample <- "TB40_siCTRL_1"
TB40CTRL_72h_2_residue$sample <- "TB40_siCTRL_2"
KSHV_72h_1_residue$sample <- "KSHV_iSLK"
NHDFRN7SK_residue$sample <- "NHDFRN7SK"
IVTRN7SK_residue$sample <- "IVTRN7SK"
HSV2_10h_residue$sample <- "HSV2"
Dumas_96h_residue$sample <- "VZV"

# Merge
merged_data_reclass_all <- rbind(
  merged_TB40_reclass_24h,
  merged_TB40_reclass_48h,
  merged_TB40_reclass_72h,
  merged_KSHV_reclass_72h_1,
  merged_NHDFRN7SK_reclass,
  merged_IVTRN7SK_reclass,
  merged_TB40CTRL_reclass_72h_1,
  merged_TB40CTRL_reclass_72h_2,
  merged_HSV2_reclass_10h,
  merged_Dumas_reclass_96h)

merged_data_reclass_summarized_all <- rbind(
  merged_TB40_reclass_summarized_24h,
  merged_TB40_reclass_summarized_48h,
  merged_TB40_reclass_summarized_72h,
  merged_KSHV_reclass_summarized_72h_1,
  merged_NHDFRN7SK_reclass_summarized,
  merged_IVTRN7SK_reclass_summarized,
  merged_TB40CTRL_reclass_summarized_72h_1,
  merged_TB40CTRL_reclass_summarized_72h_2,
  merged_HSV2_reclass_summarized_10h,
  merged_Dumas_reclass_summarized_96h)


residue_data_reclass_all <- rbind(TB40_24h_residue,
                                  TB40_48h_residue,
                                  TB40_72h_residue,
                                  TB40CTRL_72h_1_residue,
                                  TB40CTRL_72h_2_residue,
                                  KSHV_72h_1_residue,
                                  NHDFRN7SK_residue,
                                  IVTRN7SK_residue,
                                  HSV2_10h_residue,
                                  Dumas_96h_residue)


###############################################################################/



###############################################################################\
## Generating Barplots and Pie plots                                        ####
###############################################################################/

### Barplots: Decoration proportions (viral) ###################################
p_bar1 <- merged_data_reclass_all %>% 
  filter(group != "human",
         group != "IVT_RN7SK") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  mutate(group = str_to_upper(group)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("10h",
                                                  "24h",
                                                  "48h",
                                                  "72h",
                                                  "72h 1",
                                                  "72h 2",
                                                  "72h 3",
                                                  "96h"))) %>% 
  group_by(group, timepoint, class) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(group, timepoint) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ggplot(aes(x = timepoint,
             y = count,
             color = class,
             fill = class)) + 
  geom_col(position = "fill") +
  scale_color_manual(labels = c("clean", "mixed"),
                     values = c(color_bla_vir2, color_dec2)) +
  scale_fill_manual(labels = c("clean", "mixed"),
                    values = c(color_bla_vir, color_dec)) +
  geom_label(aes(label = scales::percent(percentage, accuracy = .1),
                 color = I(ifelse(class == "blank", "black", "white"))), 
             position = position_fill(vjust = .5),
             size = labelsize/.pt,
             show.legend = FALSE) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  ylab("rates of mixed tailing in herpesviral RNAs") +
  xlab("") +
  theme_clean(textbasesize) +
  my_theme_rel +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"))


### Barplots: Decoration proportions (viral) (no ncRNAs) #######################
p_bar1_no_nc <- merged_data_reclass_all %>% 
  filter(group != "human",
         group != "IVT_RN7SK") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(!str_detect(contig, "^RNA"),
         contig != "T1.1.PAN") %>% 
  mutate(group = str_to_upper(group)) %>% 
  mutate(timepoint = factor(timepoint, levels = c("10h",
                                                  "24h",
                                                  "48h",
                                                  "72h",
                                                  "72h 1",
                                                  "72h 2",
                                                  "72h 3",
                                                  "96h"))) %>% 
  group_by(group, timepoint, class) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(group, timepoint) %>% 
  mutate(percentage = count / sum(count)) %>% 
  ggplot(aes(x = timepoint,
             y = count,
             color = class,
             fill = class)) + 
  geom_col(position = "fill") +
  scale_color_manual(labels = c("clean", "mixed"),
                     values = c(color_bla_vir2, color_dec2)) +
  scale_fill_manual(labels = c("clean", "mixed"),
                    values = c(color_bla_vir, color_dec)) +
  geom_label(aes(label = scales::percent(percentage, accuracy = .1),
                 color = I(ifelse(class == "blank", "black", "white"))), 
             position = position_fill(vjust = .5),
             size = labelsize/.pt,
             show.legend = FALSE) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  ylab("rates of mixed tailing in herpesviral mRNAs") +
  xlab("") +
  theme_clean(textbasesize) +
  my_theme_rel +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"))


### Barplots: Decoration proportions (host) ####################################
p_bar2 <- merged_data_reclass_all %>% 
  filter(group != "hcmv",
         group != "kshv",
         group != "HSV2",
         group != "VZV") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  mutate(timepoint = factor(timepoint, levels = c("Mock",
                                                  "10h",
                                                  "24h",
                                                  "48h",
                                                  "72h",
                                                  "72h 1",
                                                  "72h 2",
                                                  "72h 3",
                                                  "96h"))) %>% 
  group_by(group, timepoint, class) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  group_by(group, timepoint) %>% 
  mutate(percentage = count / sum(count)) %>% 
  mutate(group = case_when(
    group == "IVT_RN7SK" ~ "IVT RN7SK",
    TRUE ~ group
  )) %>% 
  ggplot(aes(x = timepoint,
             y = count,
             color = class,
             fill = class)) + 
  geom_col(position="fill") +
  scale_color_manual(labels = c("clean", "mixed"),
                     values = c(color_bla_hum2, color_dec2)) +
  scale_fill_manual(labels = c("clean", "mixed"),
                    values = c(color_bla_hum, color_dec)) +
  geom_label(aes(label = scales::percent(percentage, accuracy = .1),
                 color = I(ifelse(class == "blank", "black", "white"))), 
             position = position_fill(vjust = .5),
             size = labelsize/.pt,
             show.legend = FALSE) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  scale_x_discrete(labels = function(x) ifelse(is.na(x), "", x)) +
  ylab("rates of mixed tailing in human mRNAs") +
  xlab("") +
  theme_clean(textbasesize) +
  my_theme_rel +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"))


### Pie plots: Decoration types - HCMV #########################################
p_pie1 <- merged_data_reclass_all %>% 
  filter(class == "decorated",
         group == "hcmv") %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  group_by(timepoint) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = row_number() <= 3) %>% 
  ungroup() %>% 
  mutate(timepoint = factor(timepoint,
                            levels = c("", "24h", "48h", "72h 1"))) %>%
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("HCMV") +
  labs(fill = "Incorporations") +
  facet_wrap(~ timepoint2,
             ncol = 2,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


### Pie plots: Decoration types - Human ########################################
p_pie2 <- merged_data_reclass_all %>% 
  filter(class == "decorated",
         group == "human") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(sample != "TB40_siCTRL_1",
         sample != "TB40_siCTRL_2") %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  group_by(timepoint) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = row_number() <= 3) %>% 
  ungroup() %>% 
  mutate(timepoint = factor(timepoint, levels = c("Mock",
                                                  "24h",
                                                  "48h",
                                                  "72h 1"))) %>%
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("human") +
  labs(fill = "Incorporations") +
  facet_wrap(~ timepoint) +
  theme_void(textbasesize) +
  my_theme_rel_pie
theme(aspect.ratio = 1)


### Pie plots: Decoration types - Other Viruses ################################
p_pie3 <- merged_data_reclass_all %>% 
  filter(class == "decorated",
         group == "kshv" |
           group == "HSV2" |
           group == "VZV") %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  mutate(group2 = factor(paste0(case_when(group == "kshv" ~ "KSHV",
                                          TRUE ~ group),
                                " (", timepoint, ")"),
                         levels = c("HSV2 (10h)",
                                    "KSHV (72h)",
                                    "VZV (96h)"))) %>% 
  group_by(group2) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = row_number() <= 3) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("other herpesviruses") +
  labs(fill = "Incorporations") +
  facet_wrap(~ group2,
             ncol = 1,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


### Pie plots: Decoration types - All Viruses ##################################
p_pie4 <- merged_data_reclass_all %>% 
  filter(class == "decorated",
         group == "hcmv" |
           group == "kshv" |
           group == "HSV2" |
           group == "VZV") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(sample != "TB40_siCTRL_1",
         sample != "TB40_siCTRL_2") %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  mutate(group2 = factor(paste0(case_when(group == "kshv" ~ "KSHV",
                                          group == "hcmv" ~ "HCMV",
                                          TRUE ~ group),
                                " (", timepoint, ")"),
                         levels = c("HCMV (24h)",
                                    "HCMV (48h)",
                                    "HCMV (72h 1)",
                                    "HSV2 (10h)",
                                    "KSHV (72h)",
                                    "VZV (96h)"))) %>% 
  group_by(group2) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = (label == "C" | label == "G" | label == "U")) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("herpesviruses") +
  labs(fill = "Incorporations") +
  facet_wrap(~ group2,
             ncol = 3,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


### Pie plots: Decoration types - All Viruses (no ncRNAs) ######################
p_pie4_no_nc <- merged_data_reclass_all %>% 
  filter(class == "decorated",
         group == "hcmv" |
           group == "kshv" |
           group == "HSV2" |
           group == "VZV") %>% 
  filter(!str_detect(contig, "^RNA"),
         contig != "T1.1.PAN") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(sample != "TB40_siCTRL_1",
         sample != "TB40_siCTRL_2") %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  mutate(group2 = factor(paste0(case_when(group == "kshv" ~ "KSHV",
                                          group == "hcmv" ~ "HCMV",
                                          TRUE ~ group),
                                " (", timepoint, ")"),
                         levels = c("HCMV (24h)",
                                    "HCMV (48h)",
                                    "HCMV (72h 1)",
                                    "HSV2 (10h)",
                                    "KSHV (72h)",
                                    "VZV (96h)"))) %>% 
  group_by(group2) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = (label == "C" | label == "G" | label == "U")) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("herpesviruses (no ncRNAs)") +
  labs(fill = "Incorporations") +
  facet_wrap(~ group2,
             ncol = 3,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


### Pie plots: Decoration types - HCMV siCTRL 72h & IVT RN7SK ##################
p_pie5 <- merged_data_reclass_all %>% 
  filter(sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2" |
           sample == "IVTRN7SK",
         class == "decorated",
         group == "hcmv" |
           group == "human" |
           group == "IVT_RN7SK") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  mutate(group2 = factor(case_when(timepoint == "72h 2" &
                                     group == "human" ~ "human (72h 2)",
                                   timepoint == "72h 3" &
                                     group == "human" ~ "human (72h 3)",
                                   timepoint == "72h 2" &
                                     group == "hcmv" ~ "HCMV (72h 2)",
                                   timepoint == "72h 3" &
                                     group == "hcmv" ~ "HCMV (72h 3)",
                                   sample == "IVTRN7SK" ~ "IVT RN7SK",
                                   TRUE ~ group),
                         levels = c("IVT RN7SK",
                                    "human (72h 2)",
                                    "human (72h 3)",
                                    "",
                                    "HCMV (72h 2)",
                                    "HCMV (72h 3)"))) %>% 
  group_by(group2) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = (label == "C" | label == "G" | label == "U")) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("HCMV 72h replicates & IVT") +
  labs(fill = "Incorporations") +
  facet_wrap(~ group2,
             ncol = 3,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


### Pie plots: Decoration types - HCMV siCTRL 72h & IVT RN7SK (no ncRNAS) ######
p_pie5_no_nc <- merged_data_reclass_all %>% 
  filter(sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2" |
           sample == "IVTRN7SK",
         class == "decorated",
         group == "hcmv" |
           group == "human" |
           group == "IVT_RN7SK") %>% 
  filter(!str_detect(contig, "^RNA"),
         contig != "T1.1.PAN") %>% 
  mutate(timepoint = case_when(
    sample == "TB40_POLYA_72" ~ "72h 1",
    sample == "TB40_siCTRL_1" ~ "72h 2",
    sample == "TB40_siCTRL_2" ~ "72h 3",
    TRUE ~ timepoint
  )) %>% 
  filter(!is.na(nonA_residues)) %>% 
  mutate(sum_vals = prediction_C + prediction_G + prediction_U) %>%
  rowwise() %>%
  mutate(label = {
    parts <- c(rep("C", prediction_C),
               rep("G", prediction_G),
               rep("U", prediction_U))
    paste(parts, collapse = " + ")
  }) %>%
  ungroup() %>%
  mutate(label = ifelse(sum_vals <= 2, label, "≥3 incorporations")) %>% 
  mutate(group2 = factor(case_when(timepoint == "72h 2" &
                                     group == "human" ~ "human (72h 2)",
                                   timepoint == "72h 3" &
                                     group == "human" ~ "human (72h 3)",
                                   timepoint == "72h 2" &
                                     group == "hcmv" ~ "HCMV (72h 2)",
                                   timepoint == "72h 3" &
                                     group == "hcmv" ~ "HCMV (72h 3)",
                                   sample == "IVTRN7SK" ~ "IVT RN7SK",
                                   TRUE ~ group),
                         levels = c("IVT RN7SK",
                                    "human (72h 2)",
                                    "human (72h 3)",
                                    "",
                                    "HCMV (72h 2)",
                                    "HCMV (72h 3)"))) %>% 
  group_by(group2) %>% 
  count(label) %>% 
  mutate(Proportion = n / sum(n)) %>%
  arrange(desc(Proportion)) %>%
  mutate(show_label = (label == "C" | label == "G" | label == "U")) %>% 
  ungroup() %>% 
  mutate(label = factor(label, levels = c("C", "U", "G", "C + C", "U + U",
                                          "G + G", "C + U", "C + G", "G + U",
                                          "≥3 incorporations"))) %>% 
  ggplot(aes(x = "", y = Proportion, fill = label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y",
              direction = -1) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_viridis_d(option = "viridis",
                       direction = -1) +  
  geom_text(aes(label = ifelse(show_label,
                               paste0(label, " (",
                                      scales::percent(Proportion,
                                                      accuracy = .1), ")"), ""),
                x = 1.2),
            position = position_stack(vjust = 0.5),
            size = labelsize/.pt) +
  ggtitle("HCMV 72h replicates & IVT (no ncRNAs)") +
  labs(fill = "Incorporations") +
  facet_wrap(~ group2,
             ncol = 3,
             drop = FALSE) +
  theme_void(textbasesize) +
  my_theme_rel_pie +
  theme(aspect.ratio = 1)


###############################################################################/



###############################################################################\
## Data rangling for Manhatten-like plots, Boxplots & Statistics            ####
###############################################################################/

### Processing for Manhatten-like, Boxplots & Statistics  ######################

## Get C/G/U frequencies
# (nonA-events per total sum of poly(A)-tail length in each
# transcript (per timepoint))
merged_data_reclass_summarized_all$C_freq <-
  merged_data_reclass_summarized_all$hits_C /
  (merged_data_reclass_summarized_all$polya_mean *
     merged_data_reclass_summarized_all$counts_total)
merged_data_reclass_summarized_all$G_freq <-
  merged_data_reclass_summarized_all$hits_G /
  (merged_data_reclass_summarized_all$polya_mean *
     merged_data_reclass_summarized_all$counts_total)
merged_data_reclass_summarized_all$U_freq <-
  merged_data_reclass_summarized_all$hits_U /
  (merged_data_reclass_summarized_all$polya_mean *
     merged_data_reclass_summarized_all$counts_total)

# Get 50 mostly expressed human Mock transcripts
most_50_Mock <- merged_data_reclass_all %>% 
  filter(timepoint == "Mock",
         group == "human") %>% 
  group_by(ensembl_transcript_id_full) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  filter(count >= 30) %>% # doesn't technically take effect,
  # as top 50 transcripts have higher depth anyway
  arrange(desc(count)) %>% 
  filter(row_number() <= 50) %>% 
  select(ensembl_transcript_id_full) %>% 
  unlist() %>% 
  unique()

# Get most expressed HCMV transcripts
# (sum of all 50 most expressed transcripts across all HCMV datasets)
most_HCMV <- merged_data_reclass_all %>% 
  filter(group == "hcmv") %>% 
  group_by(sample, ensembl_transcript_id_full) %>% 
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  filter(row_number() <= 50) %>% 
  ungroup() %>% 
  filter(count >= 30) %>% # doesn't technically take effect,
  # as top 50 transcripts have higher depth anyway
  select(ensembl_transcript_id_full) %>% 
  unlist() %>% 
  unique()


# Extract only transcripts above previously determined threshold
merged_data_reclass_summarized_all_filt <- merged_data_reclass_summarized_all %>%
  filter(group != "IVT_RN7SK",
         group != "kshv",
         group != "HSV2",
         group != "VZV") %>%
  filter(ensembl_transcript_id_full %in% most_HCMV |
           ensembl_transcript_id_full %in% most_50_Mock) %>%
  mutate(cond_col = interaction(group, timepoint, sep = ": ")) %>% 
  mutate(timepoint = factor(timepoint, levels = c("Mock",
                                                  "24h",
                                                  "48h",
                                                  "72h"))) %>%
  mutate(cond_col = factor(cond_col, levels = c("hcmv: 24h",
                                                "hcmv: 48h",
                                                "hcmv: 72h",
                                                "human: Mock",
                                                "human: 24h",
                                                "human: 48h",
                                                "human: 72h")))


mean_C_human_mock <- merged_data_reclass_summarized_all_filt %>% 
  filter(cond_col == "human: Mock") %>%
  summarise(avg = mean(C_freq, na.rm = TRUE)) %>%
  pull(avg)
mean_G_human_mock <- merged_data_reclass_summarized_all_filt %>% 
  filter(cond_col == "human: Mock") %>%
  summarise(avg = mean(G_freq, na.rm = TRUE)) %>%
  pull(avg)
mean_U_human_mock <- merged_data_reclass_summarized_all_filt %>% 
  filter(cond_col == "human: Mock") %>%
  summarise(avg = mean(U_freq, na.rm = TRUE)) %>%
  pull(avg)



### GLM Statistics #############################################################

# filter dataframe and set thresholds
merged_data_reclass_hcmv_human <- merged_data_reclass_all %>% 
  filter(group == "human" | group == "hcmv") %>% 
  mutate(polya_round = round(polya_length))

to_test_timepopint <- "72h" # latest timepoint, this will be the oe tested
min_number_reads <- 50 # min total reads across latest replicates per transcript


## Apply GLM on mixed poly(A) tailing for C
i <- 0
results_list_C <- apply_glm(merged_data_reclass_hcmv_human,
                            to_test_timepopint,
                            prediction_C)[[1]]
results_df_C <- apply_glm(merged_data_reclass_hcmv_human,
                          to_test_timepopint,
                          prediction_C)[[2]]

saveRDS(results_df_C, "./data/pooled_glm_stat_C.rds")


## Apply GLM on mixed poly(A) tailing for G
i <- 0
results_list_G <- apply_glm(merged_data_reclass_hcmv_human,
                            to_test_timepopint,
                            prediction_G)[[1]]
results_df_G <- apply_glm(merged_data_reclass_hcmv_human,
                          to_test_timepopint,
                          prediction_G)[[2]]

saveRDS(results_df_G, "./data/pooled_glm_stat_G.rds")

## Apply GLM on mixed poly(A) tailing for U
i <- 0
results_list_U <- apply_glm(merged_data_reclass_hcmv_human,
                            to_test_timepopint,
                            prediction_U)[[1]]
results_df_U <- apply_glm(merged_data_reclass_hcmv_human,
                          to_test_timepopint,
                          prediction_U)[[2]]

saveRDS(results_df_U, "./data/pooled_glm_stat_U.rds")


# get significant results only and get transcript names

results_df_C <- readRDS("./data/pooled_glm_stat_C.rds")
results_df_G <- readRDS("./data/pooled_glm_stat_G.rds")
results_df_U <- readRDS("./data/pooled_glm_stat_U.rds")


group_map <- merged_data_reclass_hcmv_human %>% 
  distinct(ensembl_transcript_id_full, group) %>% 
  rename(transcript_id = ensembl_transcript_id_full)

# Annotate dataframes with nucleotides
results_df_C <- annot_glm_df(results_df_C, "C")
results_df_G <- annot_glm_df(results_df_G, "G")
results_df_U <- annot_glm_df(results_df_U, "U")

# Get signifant transcripts
results_df_C_sig <- get_GLM_sign(results_df_C, "C")
results_df_G_sig <- get_GLM_sign(results_df_G, "G")
results_df_U_sig <- get_GLM_sign(results_df_U, "U")

#saveRDS(results_df_C_sig, "./data/pooled_glm_stat_C_sig.rds")
#saveRDS(results_df_G_sig, "./data/pooled_glm_stat_G_sig.rds")
#saveRDS(results_df_U_sig, "./data/pooled_glm_stat_U_sig.rds")


# Merge dataframes to one
results_df_all <- rbind(results_df_C , results_df_G , results_df_U)

results_df_all_sig <- rbind(distinct(results_df_C_sig, transcript_id, q_value, dec_type),
                            distinct(results_df_G_sig, transcript_id, q_value, dec_type),
                            distinct(results_df_U_sig, transcript_id, q_value, dec_type)) %>%
  group_by(transcript_id) %>% 
  spread(dec_type, q_value)


library(tidyr)

results_df_all_sig <- bind_rows(results_df_C_sig,
                                results_df_G_sig,
                                results_df_U_sig) %>%
  select(transcript_id, dec_type, q_value) %>%
  pivot_wider(names_from = dec_type, values_from = q_value)


###############################################################################/



###############################################################################\
## Generating Boxplots & Lineplots                                          ####
###############################################################################/

# mind number of reads per transcripts
bxplt_min_reads <- 30

### Boxplot (viral) ############################################################

bxplt_x_labs <- c("24h", "48h", "72h 1", "72h 2", "72h 3")

merged_data_reclass_summarized_all_bxplt <- merged_data_reclass_summarized_all %>% 
  filter(group == "hcmv",
         counts_total >= bxplt_min_reads) %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  ))

# Check for differences via Kruskall-Wallis + Dunn with Bejamini Hochberg
KW_test_C <- kruskal.test(C_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt)
KW_test_G <- kruskal.test(G_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt)
KW_test_U <- kruskal.test(U_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt)

Dunn_test_C <- dunnTest(C_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt, method="bonferroni")$res
Dunn_test_G <- dunnTest(G_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt, method="bonferroni")$res
Dunn_test_U <- dunnTest(U_freq ~ sample, data = merged_data_reclass_summarized_all_bxplt, method="bonferroni")$res

Dunn_sig_C <- subset(Dunn_test_C, P.adj <= 0.05)
Dunn_sig_G <- subset(Dunn_test_G, P.adj <= 0.05)
Dunn_sig_U <- subset(Dunn_test_U, P.adj <= 0.05)
comparisons_C <- strsplit(Dunn_sig_C$Comparison, " - ")
comparisons_G <- strsplit(Dunn_sig_G$Comparison, " - ")
comparisons_U <- strsplit(Dunn_sig_U$Comparison, " - ")

## C
set.seed(42)
p_box1 <- merged_data_reclass_summarized_all_bxplt %>%
  ggplot(aes(x = sample,
             y = 1000*(C_freq),
             color = sample,
             fill = sample)) +
  geom_jitter(size = boxplot_pointsize/.pt,
              position = position_jitterdodge(0.7),
              show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               color = "black",
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_C,
              annotations = paste0(sprintf("|Z|≈%.1f", abs(Dunn_sig_C$Z))
                                   #, " (", sprintf("p≈%.1e", Dunn_sig_C$P.adj), ")"
              ),
              y_position = seq(1.5, 2.1, length.out = nrow(Dunn_sig_C)),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt) +
  xlab("") +
  ylab("number of Cs\nper 1000 nucleotides") +
  scale_color_manual(values = c("TB40_POLYA_24" = color_vir_24,
                                "TB40_POLYA_48" = color_vir_48,
                                "TB40_POLYA_72" = color_vir_72,
                                "TB40_siCTRL_1" = color_vir_72_si,
                                "TB40_siCTRL_2" = color_vir_72_si)) +
  scale_x_discrete(labels = bxplt_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## G
set.seed(42)
p_box2 <- merged_data_reclass_summarized_all_bxplt %>%
  ggplot(aes(x = sample,
             y = 1000*(G_freq),
             color = sample,
             fill = sample)) +
  geom_jitter(size = boxplot_pointsize/.pt,
              position = position_jitterdodge(0.7),
              show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               color = "black",
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_G,
              annotations = paste0(sprintf("|Z|≈%.1f", abs(Dunn_sig_G$Z))
                                   #, " (", sprintf("p≈%.1e", Dunn_sig_G$P.adj), ")"
              ),
              y_position = seq(1.9, 2.1, length.out = nrow(Dunn_sig_G)),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt) +
  xlab("") +
  ylab("number of Gs\nper 1000 nucleotides") +
  scale_color_manual(values = c("TB40_POLYA_24" = color_vir_24,
                                "TB40_POLYA_48" = color_vir_48,
                                "TB40_POLYA_72" = color_vir_72,
                                "TB40_siCTRL_1" = color_vir_72_si,
                                "TB40_siCTRL_2" = color_vir_72_si)) +
  scale_x_discrete(labels= bxplt_x_labs) +
  coord_cartesian(ylim = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## U
set.seed(42)
p_box3 <- merged_data_reclass_summarized_all_bxplt %>%
  ggplot(aes(x = sample,
             y = 1000*(U_freq),
             color = sample,
             #color = hits_U,
             fill = sample)) +
  geom_jitter(size = boxplot_pointsize/.pt,
              position = position_jitterdodge(0.7),
              show.legend = FALSE) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               color = "black",
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_U,
              annotations = paste0(sprintf("|Z|≈%.1f", abs(Dunn_sig_U$Z))
                                   #, " (", sprintf("p≈%.1e", Dunn_sig_U$P.adj), ")"
              ),
              y_position = seq(1.3, 2.1, length.out = nrow(Dunn_sig_U)),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt) +
  xlab("") +
  ylab("number of Us\nper 1000 nucleotides") +
  scale_color_manual(values = c("TB40_POLYA_24" = color_vir_24,
                                "TB40_POLYA_48" = color_vir_48,
                                "TB40_POLYA_72" = color_vir_72,
                                "TB40_siCTRL_1" = color_vir_72_si,
                                "TB40_siCTRL_2" = color_vir_72_si)) +
  scale_x_discrete(labels= bxplt_x_labs) +
  coord_cartesian(ylim = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


### Boxplot (human) ############################################################

bxplt_hum_x_labs <- c("Mock", "24h", "48h", "72h 1", "72h 2", "72h 3")

merged_data_reclass_summarized_all_bxplt_hum <- merged_data_reclass_summarized_all %>% 
  filter(sample == "TB40_POLYA_24" |
           sample == "TB40_POLYA_48" |
           sample == "TB40_POLYA_72" |
           sample == "NHDFRN7SK" |
           sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2") %>% 
  filter(group == "human",
         counts_total >= bxplt_min_reads) %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  mutate(sample = factor(sample,
                         levels = c("NHDFRN7SK",
                                    "TB40_POLYA_24",
                                    "TB40_POLYA_48",
                                    "TB40_POLYA_72",
                                    "TB40_siCTRL_1",
                                    "TB40_siCTRL_2")))


# Check for differences via Kruskall-Wallis + Dunn with Bejamini Hochberg
KW_test_C_hum <- kruskal.test(C_freq ~ sample,
                              data = merged_data_reclass_summarized_all_bxplt_hum)
KW_test_G_hum <- kruskal.test(G_freq ~ sample,
                              data = merged_data_reclass_summarized_all_bxplt_hum)
KW_test_U_hum <- kruskal.test(U_freq ~ sample,
                              data = merged_data_reclass_summarized_all_bxplt_hum)

Dunn_test_C_hum <- dunnTest(C_freq ~ sample,
                            data = merged_data_reclass_summarized_all_bxplt_hum,
                            method="bonferroni")$res
Dunn_test_G_hum <- dunnTest(G_freq ~ sample,
                            data = merged_data_reclass_summarized_all_bxplt_hum,
                            method="bonferroni")$res
Dunn_test_U_hum <- dunnTest(U_freq ~ sample,
                            data = merged_data_reclass_summarized_all_bxplt_hum,
                            method="bonferroni")$res

Dunn_sig_C_hum <- subset(Dunn_test_C_hum, P.adj <= 0.05)
Dunn_sig_G_hum <- subset(Dunn_test_G_hum, P.adj <= 0.05)
Dunn_sig_U_hum <- subset(Dunn_test_U_hum, P.adj <= 0.05)
comparisons_C_hum <- strsplit(Dunn_sig_C_hum$Comparison, " - ")
comparisons_G_hum <- strsplit(Dunn_sig_G_hum$Comparison, " - ")
comparisons_U_hum <- strsplit(Dunn_sig_U_hum$Comparison, " - ")

# Order entries
# f_order_entries <- function(df) {
#   df$g1 <- sub(" - .*", "", df$Comparison)
#   df$g2 <- sub(".* - ", "", df$Comparison)
#   
#   lvl <- c("Mock", "24h", "48h", "72h 1", "72h 2", "72h 3")
#   df$g1 <- factor(df$g1, levels = lvl)
#   df$g2 <- factor(df$g2, levels = lvl)
#   
#   df_ordered <- df[order(df$g2, df$g1), ]
#   row.names(df_ordered) <- NULL
#   
#   df_ordered$Comparison <- paste(df_ordered$g1, df_ordered$g2, sep = " - ")
#   
#   return(df)
# }
# Dunn_sig_C_hum <- f_order_entries(Dunn_sig_C_hum)


box_max_y <- 2.2

## C
set.seed(42)
p_box1_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(C_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_C_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_C_hum)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_C_hum)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Cs\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## G
set.seed(42)
p_box2_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(G_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_G_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_G_hum)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_G_hum)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Gs\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## U
set.seed(42)
p_box3_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(U_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_U_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_U_hum)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_U_hum)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Us\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



### Violin (human) #############################################################

violin_max_y <- 11.5

## C
set.seed(42)
p_vio1_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(C_freq),
             color = sample,
             fill = sample)) +
  geom_violin(show.legend = FALSE,
              linewidth = boxplot_linesize/.pt,
              width = 1) +
  stat_summary(fun = mean,
             geom = "point",
             color = "green",
             size = boxplot_pointsize/.pt*1.35,
             position = position_dodge(1),
             na.rm = TRUE,
             show.legend = FALSE) +
  geom_signif(comparisons = comparisons_C_hum,
              annotations = "",
              y_position = seq(0.95*log(violin_max_y, base = 2) -
                                 0.025*log(violin_max_y, base = 2) *
                                 (nrow(Dunn_sig_C_hum)-1),
                               0.95*log(violin_max_y, base = 2),
                               length.out = nrow(Dunn_sig_C_hum)),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Cs\nper 1000 nucleotides") +
  scale_color_manual(values = rep("black", 6)) +
  scale_fill_manual(values = c("NHDFRN7SK" = color_hum_M,
                               "TB40_POLYA_24" = color_hum_24,
                               "TB40_POLYA_48" = color_hum_48,
                               "TB40_POLYA_72" = color_hum_72,
                               "TB40_siCTRL_1" = color_hum_72_si,
                               "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, violin_max_y),
                     transform = scales::pseudo_log_trans(base = 2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## G
set.seed(42)
p_vio2_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(G_freq),
             color = sample,
             fill = sample)) +
  geom_violin(show.legend = FALSE,
              linewidth = boxplot_linesize/.pt,
              width = 1) +
  stat_summary(fun = mean,
             geom = "point",
             color = "green",
             size = boxplot_pointsize/.pt*1.35,
             position = position_dodge(1),
             na.rm = TRUE,
             show.legend = FALSE) +
  geom_signif(comparisons = comparisons_G_hum,
              annotations = "",
              y_position = seq(0.95*log(violin_max_y, base = 2) -
                                 0.025*log(violin_max_y, base = 2) *
                                 (nrow(Dunn_sig_G_hum)-1),
                               0.95*log(violin_max_y, base = 2),
                               length.out = nrow(Dunn_sig_G_hum)),
              #y_position = seq(0.5, 1, length.out = 6),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Gs\nper 1000 nucleotides") +
  scale_color_manual(values = rep("black", 6)) +
  scale_fill_manual(values = c("NHDFRN7SK" = color_hum_M,
                               "TB40_POLYA_24" = color_hum_24,
                               "TB40_POLYA_48" = color_hum_48,
                               "TB40_POLYA_72" = color_hum_72,
                               "TB40_siCTRL_1" = color_hum_72_si,
                               "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(
    limits = c(0, 11.5),
    transform = scales::pseudo_log_trans(base = 2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## U
set.seed(42)
p_vio3_hum <- merged_data_reclass_summarized_all_bxplt_hum %>%
  ggplot(aes(x = sample,
             y = 1000*(U_freq),
             color = sample,
             fill = sample)) +
  geom_violin(show.legend = FALSE,
              linewidth = boxplot_linesize/.pt,
              width = 1) +
  stat_summary(fun = mean,
             geom = "point",
             color = "green",
             size = boxplot_pointsize/.pt*1.35,
             position = position_dodge(1),
             na.rm = TRUE,
             show.legend = FALSE) +
  geom_signif(comparisons = comparisons_U_hum,
              annotations = "",
              y_position = seq(0.95*log(violin_max_y, base = 2) -
                                 0.025*log(violin_max_y, base = 2) *
                                 (nrow(Dunn_sig_U_hum)-1),
                               0.95*log(violin_max_y, base = 2),
                               length.out = nrow(Dunn_sig_U_hum)),
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Us\nper 1000 nucleotides") +
  scale_color_manual(values = rep("black", 6)) +
  scale_fill_manual(values = c("NHDFRN7SK" = color_hum_M,
                               "TB40_POLYA_24" = color_hum_24,
                               "TB40_POLYA_48" = color_hum_48,
                               "TB40_POLYA_72" = color_hum_72,
                               "TB40_siCTRL_1" = color_hum_72_si,
                               "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, violin_max_y),
                     transform = scales::pseudo_log_trans(base = 2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


### Boxplot (human; only shared transcripts) ###################################

box_n_groups <- n_distinct(merged_data_reclass_summarized_all_bxplt_hum$sample)
merged_data_reclass_summarized_all_bxplt_hum_shared <- merged_data_reclass_summarized_all_bxplt_hum %>% 
  group_by(ensembl_transcript_id_full) %>%
  filter(n_distinct(sample) == box_n_groups) %>%
  ungroup()

# Check for differences via Kruskall-Wallis + Dunn with Bejamini Hochberg
KW_test_C_hum_shared <- kruskal.test(C_freq ~ sample,
                                     data = merged_data_reclass_summarized_all_bxplt_hum_shared)
KW_test_G_hum_shared <- kruskal.test(G_freq ~ sample,
                                     data = merged_data_reclass_summarized_all_bxplt_hum_shared)
KW_test_U_hum_shared <- kruskal.test(U_freq ~ sample,
                                     data = merged_data_reclass_summarized_all_bxplt_hum_shared)

Dunn_test_C_hum_shared <- dunnTest(C_freq ~ sample,
                                   data = merged_data_reclass_summarized_all_bxplt_hum_shared,
                                   method="bonferroni")$res
Dunn_test_G_hum_shared <- dunnTest(G_freq ~ sample,
                                   data = merged_data_reclass_summarized_all_bxplt_hum_shared,
                                   method="bonferroni")$res
Dunn_test_U_hum_shared <- dunnTest(U_freq ~ sample,
                                   data = merged_data_reclass_summarized_all_bxplt_hum_shared,
                                   method="bonferroni")$res

Dunn_sig_C_hum_shared <- subset(Dunn_test_C_hum_shared, P.adj <= 0.05)
Dunn_sig_G_hum_shared <- subset(Dunn_test_G_hum_shared, P.adj <= 0.05)
Dunn_sig_U_hum_shared <- subset(Dunn_test_U_hum_shared, P.adj <= 0.05)
comparisons_C_hum_shared <- strsplit(Dunn_sig_C_hum_shared$Comparison, " - ")
comparisons_G_hum_shared <- strsplit(Dunn_sig_G_hum_shared$Comparison, " - ")
comparisons_U_hum_shared <- strsplit(Dunn_sig_U_hum_shared$Comparison, " - ")


## C
set.seed(42)
p_box1_hum_shared <- merged_data_reclass_summarized_all_bxplt_hum_shared %>%
  ggplot(aes(x = sample,
             y = 1000*(C_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               # color = "black",
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_C_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_C_hum_shared)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_C_hum_shared)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Cs\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## G
set.seed(42)
p_box2_hum_shared <- merged_data_reclass_summarized_all_bxplt_hum_shared %>%
  ggplot(aes(x = sample,
             y = 1000*(G_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               # color = "black",
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_G_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_G_hum_shared)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_G_hum_shared)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Gs\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## U
set.seed(42)
p_box3_hum_shared <- merged_data_reclass_summarized_all_bxplt_hum_shared %>%
  ggplot(aes(x = sample,
             y = 1000*(U_freq),
             color = sample,
             fill = sample)) +
  geom_boxplot(show.legend = FALSE,
               alpha = 0.9,
               outliers = FALSE,
               fill = NA,
               linewidth = boxplot_linesize/.pt) +
  stat_summary(fun = mean,
               geom = "point",
               color = "green",
               size = boxplot_pointsize/.pt*1.35,
               position = position_dodge(1),
               na.rm = TRUE,
               show.legend = FALSE) +
  geom_signif(comparisons = comparisons_U_hum,
              annotations = "",
              y_position = seq(box_max_y - 0.025*box_max_y*(nrow(Dunn_sig_U_hum_shared)-1),
                               box_max_y,
                               length.out = nrow(Dunn_sig_U_hum_shared)),
              step_increase = 0,
              margin_top = 0,
              tip_length = 0,
              color = "black",
              textsize = boxplot_sigsize/.pt/1.5) +
  xlab("") +
  ylab("number of Us\nper 1000 nucleotides") +
  scale_color_manual(values = c("NHDFRN7SK" = color_hum_M,
                                "TB40_POLYA_24" = color_hum_24,
                                "TB40_POLYA_48" = color_hum_48,
                                "TB40_POLYA_72" = color_hum_72,
                                "TB40_siCTRL_1" = color_hum_72_si,
                                "TB40_siCTRL_2" = color_hum_72_si)) +
  scale_x_discrete(labels = bxplt_hum_x_labs) +
  scale_y_continuous(limits = c(0, 2.2)) +
  my_theme_rel +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


### Lineplots: nonA position ###################################################


## Data wrangling
# Get to-analyse transcripts (by sample/timepoint) and count how often a certain
# nucleotide occurs for RNAs with every given poly(A)-tail length (with a
# specified max length to analyse)
# Also count commulutative occureence, which is important to later weight the
# occurrence by poly(A)-tail length (e.g. to calculate tie frequency of a 'G'
# 50 nucleotides away of the poly(A) tail 3' end, the number of Gs in the sub-
# set of RNAs wit a poly(A)-tail of >= 50 is relevant, as only these could on
# principle contain a G at 50 nt distance).
n_reads_by_polyalen <- merged_data_reclass_all %>% 
  filter(sample == "TB40_POLYA_24" |
           sample == "TB40_POLYA_48" |
           sample == "TB40_POLYA_72" |
           sample == "NHDFRN7SK" |
           sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2") %>% 
  filter(group == "human" |
           group == "hcmv") %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  mutate(group2 = factor(case_when(
    group == "human" ~ "human",
    transcript_name == "RNA2.7-1" ~ "RNA2.7",
    transcript_name == "mRNA.UL69-1" ~ "UL69",
    group == "hcmv" ~ "HCMV*",
    TRUE ~ NA
  ),
  levels = c("human", "HCMV*", "RNA2.7", "UL69"))) %>% 
  mutate(polya_length = round(polya_length)) %>% 
  group_by(timepoint, group2, polya_length) %>% 
  summarise(count = n(), .groups = "drop_last") %>%
  tidyr::complete(polya_length = 0:300, fill = list(count = 0)) %>% # set max
  # poly(A) tail length
  arrange(desc(polya_length), .by_group = TRUE) %>%
  mutate(count_ge = cumsum(count)) %>%
  ungroup()

# Calculate the distance to the poly(A)-tail 3' end based on the ninetails
# predicted position and the poly(A)-tail length of the given RNA
residue_data_reclass_all_pos <- residue_data_reclass_all %>% 
  filter(sample == "TB40_POLYA_24" |
           sample == "TB40_POLYA_48" |
           sample == "TB40_POLYA_72" |
           sample == "NHDFRN7SK" |
           sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2") %>% 
  filter(group == "human" |
           group == "hcmv") %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  mutate(group2 = factor(case_when(
    group == "human" ~ "human",
    transcript_name == "RNA2.7-1" ~ "RNA2.7",
    transcript_name == "mRNA.UL69-1" ~ "UL69",
    group == "hcmv" ~ "HCMV*",
    TRUE ~ NA
  ),
  levels = c("human", "HCMV*", "RNA2.7", "UL69"))) %>% 
  mutate(polya_length = round(polya_length),
         est_nonA_pos = round(est_nonA_pos),
         dist_3 = polya_length - est_nonA_pos)

# count nonAs
n_nonA_dist <- residue_data_reclass_all_pos %>% 
  count(timepoint, group2, dist_3, prediction, name = "n_nonA") %>%
  tidyr::complete(timepoint, group2, dist_3, prediction, fill = list(n_nonA = 0))

# Width of nucleotides to calculate sliding window frequency for mixed tailing
window_width <- 5

# Calculate nonA frequencies using the other dataframes,
# use sliding window and mean
frequency_nonA_dist <- n_nonA_dist %>%
  left_join(n_reads_by_polyalen %>% select(timepoint, group2, polya_length , count_ge),
            by = c("timepoint", "group2", "dist_3" = "polya_length")) %>% 
  group_by(timepoint, group2, prediction) %>%
  mutate(count_ge = if_else(
    is.na(count_ge) & `dist_3` < 10,
    first(count_ge[`dist_3` == 10]),
    count_ge
  )) %>%
  ungroup() %>% 
  mutate(frequency = n_nonA/count_ge) %>% 
  mutate(timepoint = factor(timepoint, levels = c("Mock", "24h", "48h", "72h"))) %>% 
  tidyr::complete(fill = list(frequency = 0)) %>% 
  filter(!(timepoint == "Mock" & group2 != "human")) %>%
  group_by(timepoint, group2, prediction) %>% 
  arrange(dist_3) %>%
  mutate(frequency_roll = 1000*rollapply(frequency, width = window_width, FUN = mean,
                                         align = "center", fill = NA, na.rm = FALSE))


# Plot
p_line <- frequency_nonA_dist %>% 
  filter(!dist_3 > 150) %>% 
  filter(!(group2 == "UL69" & (timepoint == "24h" | timepoint == "48h"))) %>% 
  ggplot(aes(x = dist_3, y = frequency_roll, colour = group2)) +
  geom_line(linewidth = 1.5*lineplot_linesize) +
  scale_x_reverse(lim = c(150, 0)) +
  facet_grid(prediction ~ timepoint,
             scales = "free_y") +
  ylab("position specific number of non-As per 1000 nucleotides") +
  xlab("distance to poly(A) 3' end") +
  scale_color_manual(values = c("human" = color_hum_M,
                              "HCMV*" = color_vir_72,
                              "RNA2.7" = color_2_7,
                              "UL69" = color_UL69)) +
  labs(color = "") +
  my_theme_rel +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))



### Lineplots: nonA position (further transcripts, collapsed) ##################

## Data wrangling
# Get to-analyse transcripts (by sample/timepoint) and count how often a certain
# nucleotide occurs for RNAs with every given poly(A)-tail length (with a
# specified max length to analyse)
# Also count commulutative occureence, which is important to later weight the
# occurrence by poly(A)-tail length (e.g. to calculate tie frequency of a 'G'
# 50 nucleotides away of the poly(A) tail 3' end, the number of Gs in the sub-
# set of RNAs wit a poly(A)-tail of >= 50 is relevant, as only these could on
# principle contain a G at 50 nt distance).
n_reads_by_polyalen_extend_coll <- merged_data_reclass_all %>% 
  filter(sample == "TB40_POLYA_24" |
           sample == "TB40_POLYA_48" |
           sample == "TB40_POLYA_72" |
           sample == "NHDFRN7SK" |
           sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2") %>% 
  filter(group == "human" |
           group == "hcmv") %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  filter(transcript_name == "RNA2.7-1" |
           transcript_name == "mRNA.UL69-1" |
           transcript_name == "mRNA.UL77-1" |
           transcript_name == "ENST00000648006.3" |
           transcript_name == "ENST00000216281.13" |
           transcript_name == "mRNA.UL22A" |
           transcript_name == "mRNA.UL132-1"
  ) %>% 
  mutate(polya_length = round(polya_length)) %>% 
  group_by(timepoint, transcript_name, polya_length) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  tidyr::complete(polya_length = 0:300, fill = list(count = 0)) %>% # set max
  # poly(A) tail length
  arrange(desc(polya_length), .by_group = TRUE) %>%
  mutate(count_ge = cumsum(count)) %>%
  ungroup()

# Calculate the distance to the poly(A)-tail 3' end based on the ninetails
# predicted position and the poly(A)-tail length of the given RNA
residue_data_reclass_all_pos_extend_coll <- residue_data_reclass_all %>% 
  filter(sample == "TB40_POLYA_24" |
           sample == "TB40_POLYA_48" |
           sample == "TB40_POLYA_72" |
           sample == "NHDFRN7SK" |
           sample == "TB40_siCTRL_1" |
           sample == "TB40_siCTRL_2") %>% 
  filter(group == "human" |
           group == "hcmv") %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  filter(transcript_name == "RNA2.7-1" |
           transcript_name == "mRNA.UL69-1" |
           transcript_name == "mRNA.UL77-1" |
           transcript_name == "ENST00000648006.3" |
           transcript_name == "ENST00000216281.13" |
           transcript_name == "mRNA.UL22A" |
           transcript_name == "mRNA.UL132-1"
  ) %>% 
  mutate(polya_length = round(polya_length),
         est_nonA_pos = round(est_nonA_pos),
         dist_3 = polya_length - est_nonA_pos)

# count nonAs
n_nonA_dist_extend_coll <- residue_data_reclass_all_pos_extend_coll %>% 
  count(timepoint, transcript_name, dist_3, prediction, name = "n_nonA") %>%
  tidyr::complete(timepoint, transcript_name, dist_3, prediction, fill = list(n_nonA = 0))

# Width of nucleotides to calculate sliding window frequency for mixed tailing
window_width <- 5

# Calculate nonA frequencies using the other dataframes,
# use sliding window and mean
frequency_nonA_dist_extend_coll <- n_nonA_dist_extend_coll %>%
  left_join(n_reads_by_polyalen_extend_coll %>% select(timepoint, transcript_name, polya_length , count_ge),
            by = c("timepoint", "transcript_name", "dist_3" = "polya_length")) %>% 
  group_by(timepoint, transcript_name, prediction) %>%
  mutate(count_ge = if_else(
    is.na(count_ge) & `dist_3` < 10,
    first(count_ge[`dist_3` == 10]),
    count_ge
  )) %>%
  ungroup() %>% 
  mutate(frequency = n_nonA/count_ge,
         group3 = interaction(timepoint, transcript_name, sep = ": ")) %>% 
  tidyr::complete(fill = list(frequency = 0)) %>% 
  group_by(timepoint, transcript_name, prediction) %>% 
  arrange(dist_3) %>%
  mutate(frequency_roll = 1000*rollapply(frequency, width = window_width, FUN = mean,
                                         align = "center", fill = NA, na.rm = FALSE))


# Plot
p_line_extend_coll <- frequency_nonA_dist_extend_coll %>% 
  filter(!dist_3 > 150) %>% 
  mutate(timepoint = factor(timepoint,
                            levels = c("Mock", "24h", "48h", "72h"))) %>% 
  filter(
    transcript_name == "ENST00000648006.3" |
      transcript_name == "ENST00000216281.13" |
      transcript_name == "mRNA.UL22A" |
      transcript_name == "mRNA.UL77-1" |
      transcript_name == "mRNA.UL132-1") %>%
  mutate(transcript_name = factor(case_when(
    transcript_name == "mRNA.UL22A" ~ "UL22A",
    transcript_name == "mRNA.UL132-1" ~ "UL132",
    transcript_name == "mRNA.UL69-1" ~ "UL69",
    transcript_name == "mRNA.UL77-1" ~ "UL77",
    transcript_name == "ENST00000648006.3" ~ "B2M-211",
    transcript_name == "ENST00000216281.13" ~ "HSP90AA1-201",
    TRUE ~ transcript_name
  ),
  level = c("UL22A", "UL132", "UL69", "UL77", "B2M-211", "HSP90AA1-201"))) %>% 
  ggplot(aes(x = dist_3, y = frequency_roll, colour = transcript_name)) +
  geom_line(linewidth = lineplot_linesize) +
  scale_x_reverse(lim = c(150, 0)) +
  facet_grid(prediction ~ timepoint,
             scales = "free_y") +
  ylab("position specific number of\nnon-As per 1000 nucleotides") +
  xlab("distance to poly(A) 3' end") +
  scale_color_manual(values = c("UL22A" = color_UL22A,
                              "UL132" = color_UL132,
                              "UL77" = color_UL77,
                              "B2M-211" = color_B2M,
                              "HSP90AA1-201" = color_HSP)) +
  labs(color = "") +
  my_theme_rel +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))


###############################################################################/



###############################################################################\
## Generating Manhatten-like plots                                          ####
###############################################################################/


### Manhatten-like-plots: nonA frequencies per transcript (hcmv + human) #######

# Create annotation
merged_data_reclass_summarized_all_filt_sig <- merged_data_reclass_summarized_all_filt %>% 
  filter(ensembl_transcript_id_full %in% most_HCMV |
           ensembl_transcript_id_full %in% most_50_Mock) %>% 
  mutate(transcript_name = case_when(
    group == "human" ~ ensembl_transcript_id_full,
    group == "hcmv" ~ gsub("(^[m]?RNA.*)::.*", "\\1", ensembl_transcript_id_full),
  )) %>% 
  rename(transcript_id = ensembl_transcript_id_full) %>% 
  left_join(results_df_all_sig, by = "transcript_id") %>% 
  mutate(
    C = case_when(
      timepoint != "72h" ~ NA,
      TRUE ~ C
    ),
    G = case_when(
      timepoint != "72h" ~ NA,
      TRUE ~ G
    ),
    U = case_when(
      timepoint != "72h" ~ NA,
      TRUE ~ U
    )) %>% 
  arrange(desc(counts_total)) %>%
  mutate(transcript_name = fct_inorder(transcript_name))

label_list <- list(
  "hcmv" = "HCMV",
  "human" = "human"
)

G_transcript_labels <- data.frame(x = c("RNA2.7-1", "mRNA.UL69-1"),
                                  lab = c("RNA2.7", "UL69"),
                                  group = "hcmv")

## C
set.seed(42)
p_man1 <- merged_data_reclass_summarized_all_filt_sig %>% 
  slice_sample(prop = 1) %>% 
  ggplot(aes(x = transcript_name,
             y = 1000*C_freq,
             color = cond_col,
             label = C
  )) +
  geom_point(size = manh_pointsize/.pt) +
  geom_text(aes(label = ifelse(timepoint == "72h" & !is.na(C), "*", "")),
            color = "black",
            y = 2.25,
            size = manh_sigsize/.pt) +
  coord_cartesian(ylim = c(0, 2.2)) +
  geom_hline(yintercept = 1000*mean_C_human_mock,
             color = "red",
             linetype = "dashed",
             linewidth = manh_linesize/.pt) +
  guides(size = "none",
         alpha = "none",
         color = guide_legend(override.aes = list(size = manh_pointsize/.pt),
                              title = "Sample")) +
  scale_x_discrete(expand = expansion(add = 1.75)) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x",
             labeller = change_facet_labels) +
  scale_color_manual(labels = c("HCMV: 24h",
                                "HCMV: 48h",
                                "HCMV: 72h",
                                "human: Mock",
                                "human: 24h",
                                "human: 48h",
                                "human: 72h"),
                     values = c("hcmv: 24h" = color_vir_24,
                                "hcmv: 48h" = color_vir_48,
                                "hcmv: 72h" = color_vir_72,
                                "human: Mock" = color_hum_M,
                                "human: 24h" = color_hum_24,
                                "human: 48h" = color_hum_48,
                                "human: 72h" = color_hum_72)) +
  labs(color = "") +
  xlab("transcripts") +
  ylab("number of Cs per 1000 nucleotides") +
  my_theme_rel +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank())

## G
set.seed(42)
p_man2 <- merged_data_reclass_summarized_all_filt_sig %>% 
  slice_sample(prop = 1) %>% 
  ggplot(aes(x = transcript_name,
             y = 1000*G_freq,
             color = cond_col,
             label = G
  )) +
  geom_point(size = manh_pointsize/.pt) +
  geom_text(aes(label = ifelse(timepoint == "72h" & !is.na(G), "*", "")),
            color = "black",
            y = 2.25,
            size = manh_sigsize/.pt) +
  geom_text(
    data = G_transcript_labels,
    inherit.aes = FALSE,
    aes(x = x, label = lab),
    y = 1.75,
    angle = 90,
    size = boxplot_sigsize/.pt) +
  coord_cartesian(ylim = c(0, 2.2)) +
  geom_hline(yintercept = 1000*mean_G_human_mock,
             color = "red",
             linetype = "dashed",
             linewidth = manh_linesize/.pt) +
  guides(size = "none",
         alpha = "none",
         color = guide_legend(override.aes = list(size = manh_pointsize/.pt),
                              title = "Sample")) +
  scale_x_discrete(expand = expansion(add = 1.75)) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x",
             labeller = change_facet_labels) +
  scale_color_manual(labels = c("HCMV: 24h",
                                "HCMV: 48h",
                                "HCMV: 72h",
                                "human: Mock",
                                "human: 24h",
                                "human: 48h",
                                "human: 72h"),
                     values = c("hcmv: 24h" = color_vir_24,
                                "hcmv: 48h" = color_vir_48,
                                "hcmv: 72h" = color_vir_72,
                                "human: Mock" = color_hum_M,
                                "human: 24h" = color_hum_24,
                                "human: 48h" = color_hum_48,
                                "human: 72h" = color_hum_72)) +
  labs(color = "") +
  xlab("transcripts") +
  ylab("number of Gs per 1000 nucleotides") +
  my_theme_rel +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank())

## U
set.seed(42)
p_man3 <- merged_data_reclass_summarized_all_filt_sig %>% 
  slice_sample(prop = 1) %>% 
  ggplot(aes(x = transcript_name,
             y = 1000*U_freq,
             color = cond_col,
             label = U
  )) +
  geom_point(size = manh_pointsize/.pt) +
  geom_text(aes(label = ifelse(timepoint == "72h" & !is.na(U), "*", "")),
            color = "black",
            y = 2.25,
            size = manh_sigsize/.pt) +
  coord_cartesian(ylim = c(0, 2.2)) +
  geom_hline(yintercept = 1000*mean_U_human_mock,
             color = "red",
             linetype = "dashed",
             linewidth = manh_linesize/.pt) +
  guides(size = "none",
         alpha = "none",
         color = guide_legend(override.aes = list(size = manh_pointsize/.pt),
                              title = "Sample")) +
  scale_x_discrete(expand = expansion(add = 1.75)) +
  facet_grid(~ group,
             scales = "free_x",
             space = "free_x",
             switch = "x",
             labeller = change_facet_labels) +
  scale_color_manual(labels = c("HCMV: 24h",
                                "HCMV: 48h",
                                "HCMV: 72h",
                                "human: Mock",
                                "human: 24h",
                                "human: 48h",
                                "human: 72h"),
                     values = c("hcmv: 24h" = color_vir_24,
                                "hcmv: 48h" = color_vir_48,
                                "hcmv: 72h" = color_vir_72,
                                "human: Mock" = color_hum_M,
                                "human: 24h" = color_hum_24,
                                "human: 48h" = color_hum_48,
                                "human: 72h" = color_hum_72)) +
  labs(color = "") +
  xlab("transcripts") +
  ylab("number of Usper 1000 nucleotides") +
  my_theme_rel +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank())


###############################################################################/



###############################################################################\
## Plotting                                                                 ####
###############################################################################/

### Collect plots ##############################################################

# extract legend to reassemble later on,
# blank some axis labels for clearer plot layouts
legend_pi <- get_legend(p_pie1)
p_pie1n <- p_pie1 + theme(legend.position = "none")
p_pie2n <- p_pie2 + theme(legend.position = "none")
p_pie3n <- p_pie3 + theme(legend.position = "none")
p_pie4n <- p_pie4 + theme(legend.position = "none")
p_pie4_no_ncn <- p_pie4_no_nc + theme(legend.position = "none")

legend_man <- get_legend(p_man1)
p_box1n <- p_box1 + theme(legend.position = "none")
p_box2n <- p_box2 + theme(legend.position = "none")
p_box3n <- p_box3 + theme(legend.position = "none")
p_man1n <- p_man1 + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_man2n <- p_man2 + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_man3n <- p_man3 + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box1n_hum <- p_box1_hum + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box2n_hum <- p_box2_hum + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box3n_hum <- p_box3_hum + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box1n_hum_shared <- p_box1_hum_shared + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box2n_hum_shared <- p_box2_hum_shared + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box3n_hum_shared <- p_box3_hum_shared + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank())
p_box1n <- p_box1n + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_box2n <- p_box2n + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_box1n_hum <- p_box1n_hum + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_box2n_hum <- p_box2n_hum + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_box1n_hum_shared <- p_box1n_hum_shared + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_box2n_hum_shared <- p_box2n_hum_shared + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_man1n <- p_man1n + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_man2n <- p_man2n + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.background = element_blank(), strip.text.x = element_blank())
p_man3n <- p_man3n + theme(axis.title.x = element_text(margin = margin(t = -(2*textbasesize))))


## Stitch plots together

# fig 6
fig6 <- (p_bar2 & labs(tag = "a")) + (p_pie2n & labs(tag = "b")) + (p_bar1 & labs(tag = "c")) + (p_pie4n & labs(tag = "d")) + legend_pi +
  plot_layout(design = "
              AAAABBBE
              CCCCDDDD
              ") &
  theme(plot.tag = element_text(size = plotlabelsize, face = "bold"))
grDevices::cairo_pdf("./plots/fig6.pdf",
                     width = 28.125,
                     height = 25)
fig6
dev.off()


# figure 7
fig7_total <- (free(side = "l", p_box1n) & labs(tag = "a")) + (p_man1n & labs(tag = "b")) + (free(side = "l", p_box1n_hum) & labs(tag = "c")) +     free(side = "l", p_box2n) + p_man2n + free(side = "l", p_box2n_hum) +     free(side = "l", p_box3n) + p_man3n + free(side = "l", p_box3n_hum) +     legend_man + (p_line & labs(tag = "d")) +
  plot_layout(design = "
              AAABBBBBBCCJJKKKKKKKKKKKKK
              DDDEEEEEEFFJJKKKKKKKKKKKKK
              GGGHHHHHHIIJJKKKKKKKKKKKKK
              ") &
  theme(plot.tag = element_text(size = plotlabelsize, face = "bold"))

grDevices::cairo_pdf("./plots/fig7.pdf",
                     width = 50,
                     height = 25)
fig7_total
dev.off()


# Suppl. figure 1
fig7_suppl2 <- p_box1_hum_shared + p_box2_hum_shared + p_box3_hum_shared +
  plot_layout(design = "
              ABC
              ")
grDevices::cairo_pdf("./plots/figS1.pdf",
                     width = 15,
                     height = 8.3)
fig7_suppl2
dev.off()


# Suppl. figure 2a/b
fig6_suppl <- p_pie5
grDevices::cairo_pdf("./plots/figS2ab.pdf",
                     width = 14,
                     height = 9)
fig6_suppl
dev.off()


# Suppl. figure 2c
fig7_suppl <- p_line_extend_coll
grDevices::cairo_pdf("./plots/figS2c.pdf",
                     width = 16.67,
                     height = 9.38)
fig7_suppl
dev.off()


###############################################################################/
