###############################################################################\
#                           Pipeline for Ninetails                          ####
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


###############################################################################/



###############################################################################\
## Loading environments ####

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(stringr)
library(dplyr)

# Set env based on your specific paths and environment names
Sys.setenv(RETICULATE_CONDA="/path/to/Anaconda3/2022.05/bin/conda")
reticulate::use_condaenv("r-ninetails")


###############################################################################/



###############################################################################\
## Run fixer script ####

# This pipeline uses Ninetails v1.0.3. To avoid an issue that can arise when
# Ninetails tries to read FAST5 files that contain more than one basecall
# groups, run the following fixer script first (and apply the changed
# functions to affected datasets).

source("./fix_multiple_basecallgroups_ninetails.R")


###############################################################################/



###############################################################################\
## Functions ####

### Run Ninetails (Basecall group 000) ####
do_ninetails_runs_0 <- function(base_path, basename, groupname, seqsum_path, workspace_path) {
  print(paste0("Now processing: ", basename, " with ", groupname))
  ninetails_run <- ninetails::check_tails(
    nanopolish = paste0(base_path, basename, ".", groupname, ".tsv"),
    sequencing_summary = paste0(seqsum_path, "sequencing_summary.txt"),
    workspace = workspace_path,
    num_cores = 8,
    basecall_group = 'Basecall_1D_000', # check for Guppys group in that fast5 files first
    pass_only = FALSE,
    save_dir = base_path)
  # Save Ninetails results in form of RDS objects as intermediate file
  saveRDS(ninetails_run, file = paste0(base_path, basename, ".", groupname, ".ninetails_pF.rds"))
  
  gc()
}


### Run Ninetails (Basecall group 001) ####
do_ninetails_runs_1 <- function(base_path, basename, groupname, seqsum_path, workspace_path) {
  print(paste0("Now processing: ", basename, " with ", groupname))
  ninetails_run <- new_check_tails(
    nanopolish = paste0(base_path, basename, ".", groupname, ".tsv"),
    sequencing_summary = paste0(seqsum_path, "sequencing_summary.txt"),
    workspace = workspace_path,
    num_cores = 8,
    basecall_group = 'Basecall_1D_001', # check for Guppys group in that fast5 files first
    pass_only = FALSE,
    save_dir = base_path)
  # Save Ninetails results in form of RDS objects as intermediate file
  saveRDS(ninetails_run, file = paste0(base_path, basename, ".", groupname, ".ninetails_pF.rds"))
  
  gc()
}


### Extract Ninetails results to usable dataframes and reclassify (increase accuracy) ####
# Reclassifies (Ninetails function) data, merges the human and viral datasets
# and writes to tsv
postproscess_ninetails <- function(base_path,
                                   basename,
                                   groupname_1,
                                   groupname_2,
                                   group_id_1,
                                   group_id_2,
                                   name_extension) {
  
  # Load RDS files
  rds_1 <- readRDS(paste0(base_path, basename, ".", groupname_1, ".ninetails_pF.rds"))
  
  # Get Ninetails class and residue tables
  class_1 <- rds_1[[1]]
  residue_1 <- rds_1[[2]]
  
  class_1$group <- group_id_1
  residue_1$group <- group_id_1
  
  class_data <- residue_data <- NULL
  
  # Apply for human subset only if it exists, otherwise work with viral only
  if(!is.null(groupname_2)){
    # Load RDS files
    rds_2 <- readRDS(paste0(base_path, basename, ".", groupname_2, ".ninetails_pF.rds"))
    
    # Get Ninetails class and residue tables
    class_2 <- rds_2[[1]]
    residue_2 <- rds_2[[2]]
    
    class_2$group <- group_id_2
    residue_2$group <- group_id_2
    
    # merge human and viral dataframes
    class_data <- rbind(class_1, class_2)
    residue_data <- rbind(residue_1, residue_2)
    
  } else {
    
    class_data <- viral_class
    residue_data <- viral_residue
    
  }
  
  # Reclassify the data
  
  # This section was taken from the Ninetails read function (excluding 'correct_labels', which seems to cause problems)
  # correct annotation <- if gencode format: create new columns with ensembl IDs
  # else created columns could be easily dropped
  # this code chunk was originally written by Paweł Krawczyk (smaegol) & incorporated in NanoTail package
  transcript_names <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", class_data$contig)
  class_data$transcript <- transcript_names
  ensembl_transcript_ids <- gsub("^(.*?)\\|.*\\|.*", "\\1", class_data$contig)
  ensembl_transcript_ids_short <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids) # without version number
  class_data$ensembl_transcript_id_full <- ensembl_transcript_ids
  class_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short
  
  transcript_names_r <- gsub(".*?\\|.*?\\|.*?\\|.*?\\|.*?\\|(.*?)\\|.*", "\\1", residue_data$contig)
  residue_data$transcript <- transcript_names_r
  ensembl_transcript_ids_r <- gsub("^(.*?)\\|.*\\|.*", "\\1", residue_data$contig)
  ensembl_transcript_ids_short_r <- gsub("(.*)\\..*", "\\1", ensembl_transcript_ids_r) # without version number
  residue_data$ensembl_transcript_id_full <- ensembl_transcript_ids_r
  residue_data$ensembl_transcript_id_short <- ensembl_transcript_ids_short_r
  
  
  ninetails_data_reclass <- ninetails::reclassify_ninetails_data(residue_data = residue_data,
                                                                 class_data = class_data,
                                                                 grouping_factor = "group", 
                                                                 transcript_column = "ensembl_transcript_id_short",
                                                                 ref = 'hsapiens')
  
  # Retrieve the data frames
  class_data_reclass <- ninetails_data_reclass[[1]]
  residue_data_reclass <- ninetails_data_reclass[[2]]
  
  # Save the data frames to tsv files
  write.table(class_data_reclass,
              file = paste0(base_path, basename, name_extension, ".class_data_reclass.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  write.table(residue_data_reclass,
              file = paste0(base_path, basename, name_extension, ".residue_data_reclass.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}




###############################################################################/



###############################################################################\
## Apply Ninetails on different datasets ####

# Note: Datasets used and how they are refernced in the R-scripts:
# HCMV (TB40/E) infected NHDF, 24 hpi, MOI = 3: TB40_24h
# HCMV (TB40/E) infected NHDF, 48 hpi, MOI = 3: TB40_48h
# HCMV (TB40/E) infected NHDF, 72 hpi 1, MOI = 3: TB40_72h
# HCMV (TB40/E) infected NHDF, 72 hpi 2 (siCTRL 1): TB40_72h_CTRL-1 / TB40CTRL_72h_1
# HCMV (TB40/E) infected NHDF, 72 hpi 3 (siCTRL 2): TB40_72h_CTRL-2 / TB40CTRL_72h_2
# Mock NHDF: NHDFpolyA-RN7SK-002 / NHDFRN7SK
#            (sample was originally created with a RN7SK spike-in, which was
#             removed downstream in silico)
# RN7SK IVT: RN7SKpolyA-002 / IVT_RN7SK
# HSV-2 infected ARPE-19, 10h hpi: 
# VZV (EMC-1) infected MeWo, 96 hpi: EMC1-MeWo-96h / EMC1-MeWo-96h-polyA / Dumas_96h
# KSHV infected iSLK, 72 hpi: KSHV-iSLK-72h / KSHV-iSLK-72h / KSHV_72h_1

# Note: abbreviations:
# tx: transcriptome
# pF: pass_only = FALSE (Ninetails parameter)
# noMT: no mitochondrial RNAs
# protEnc: protein encoding RNAs only

# Note: Ninetails runs are partially split (e.g. human and viral) for resource
# reasons. Further annotiationnecessary (e.g. timepoints) are added in the
# plotting script.

# Note: Human alignment was only included for HCMV infections.


### TB40 24h ###################################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "TB40_24h"
seqsum_path <- "../input/TB40_24h-hac.v6.1.7/"
workspace_path <- "../input/TB40_24h-hac.v6.1.7/workspace/fast5_pass/"
name_extension <- ".TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF"

# human
groupname_hum <- "GRCh38.tx.polyA.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

# hcmv
groupname_vir <- "TB40txome.polyA"
group_id_vir <- "hcmv"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### TB40 48h ###################################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "TB40_48h"
seqsum_path <- "../input/TB40_48h-hac.v6.1.7/"
workspace_path <- "../input/TB40_48h-hac.v6.1.7/workspace/20200130_1934_MN24978_FAL86847_be4e2264/fast5_all/"
name_extension <- ".TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF"

# human
groupname_hum <- "GRCh38.tx.polyA.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

# hcmv
groupname_vir <- "TB40txome.polyA"
group_id_vir <- "hcmv"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### TB40 72h ###################################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "TB40_72h"
seqsum_path <- "../input/TB40_72h-hac.v6.1.7/"
workspace_path <- "../input/TB40_72h-hac.v6.1.7/workspace/20200207_1844_MN24978_FAL81867_7f514654/fast5_all/"
name_extension <- ".TB40txome.GRCh38.tx-noMT_protEnc.ninetails_pF"

# human
groupname_hum <- "GRCh38.tx.polyA.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

# hcmv
groupname_vir <- "TB40txome.polyA"
group_id_vir <- "hcmv"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### TB40_72h_CTRL-1 ############################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "TB40_72h_CTRL-1.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf"
seqsum_path <- "../input/TB40_72h_CTRL-1-hac.v6.1.7/"
workspace_path <- "../input/TB40_72h_CTRL-1-hac.v6.1.7/workspace/20201130_1639_MN24978_FAO99310_497d8769/fast5_all/"
name_extension <- ".human-noMT_protEnc.hcmv.ninetails_pF"

# human
groupname_hum <- "polyA.human.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

# hcmv
groupname_vir <- "polyA.hcmv"
group_id_vir <- "hcmv"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### TB40_72h_CTRL-2 ############################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "TB40_72h_CTRL-2.guppy.hac.6.1.7.gc47_TB40v1.3_ENO2.tx.uf"
seqsum_path <- "../input/TB40_72h_CTRL-2-hac.v6.1.7/"
workspace_path <- "../input/TB40_72h_CTRL-2-hac.v6.1.7/workspace/20210506_1601_MN24978_FAP02780_8131a338/fast5_all/"
name_extension <- ".human-noMT_protEnc.hcmv.ninetails_pF"

# human
groupname_hum <- "polyA.human.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

# hcmv
groupname_vir <- "polyA.hcmv"
group_id_vir <- "hcmv"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### HSV2-ARPE19-10h-1 ##########################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "HSV2-ARPE19-10h-1.guppy.hac.6.1.7"
seqsum_path <- "../input//HSV2-ARPE19-10h-1-hac.v6.1.7_redo/"
workspace_path <- "../input/HSV2-ARPE19-10h-1-hac.v6.1.7_redo/workspace/20211025_1359_MN26830_FAR15445_e81919df/fast5_pass/"
name_extension <- ".HSV2-MS.tx.uf.ninetails_pF"

# human
groupname_hum <- NULL
group_id_hum <- NULL

# hsv2
groupname_vir <- "HSV2-MS.uf.polyA"
group_id_vir <- "HSV2"
do_ninetails_runs_1(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_vir,
                       groupname_2 = groupname_hum,
                       group_id_1 = group_id_vir,
                       group_id_2 = group_id_hum,
                       name_extension = name_extension)


### EMC1-MeWo-96h ##############################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "EMC1-MeWo-96h-polyA.guppy.hac.6.1.7"
seqsum_path <- "../input/EMC1-MeWo-96h-polyA-hac.v6.1.7_redo/"
workspace_path <- "../input/EMC1-MeWo-96h-polyA-hac.v6.1.7_redo/workspace/20230308_1540_MN32174_FAT75568_f6dd48b0/fast5_all/"
name_extension <- ".dumas.tx.uf.ninetails_pF"

# human
groupname_hum <- NULL
group_id_hum <- NULL

# vzv
groupname_vir <- "dumas.uf.polyA"
group_id_vir <- "VZV"
do_ninetails_runs_0(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_vir,
                       groupname_2 = groupname_hum,
                       group_id_1 = group_id_vir,
                       group_id_2 = group_id_hum,
                       name_extension = name_extension)


### KSHV-iSLK-72h-1-hac.v6.1.7_redo ############################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "KSHV-iSLK-72h-1.guppy.hac.6.1.7"
seqsum_path <- "../input/KSHV-iSLK-72h-1-hac.v6.1.7_redo/"
workspace_path <- "../input/KSHV-iSLK-72h-1-hac.v6.1.7_redo/workspace/20220919_1241_MN32174_FAT75712_ab0a3816/fast5_all/"
name_extension <- ".KSHVtxome.uf.ninetails_pF"

# human
groupname_hum <- NULL
group_id_hum <- NULL

# kshv
groupname_vir <- "KSHVtxome.uf.polyA"
group_id_vir <- "kshv"
do_ninetails_runs_0(base_path, basename, groupname, seqsum_path, workspace_path)

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_vir,
                       groupname_2 = groupname_hum,
                       group_id_1 = group_id_vir,
                       group_id_2 = group_id_hum,
                       name_extension = name_extension)


### NHDFRN7SKpolyA-002 (Mock) ##################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "NHDFpolyA-RN7SK-002.guppy.hac.6.1.7"
seqsum_path <- "../input/NHDFpolyA-RN7SK-002-hac.v6.1.7/"
workspace_path <- "../input/NHDFpolyA-RN7SK-002-hac.v6.1.7/workspace/20250908_1302_MN32174_FAT84350_7ecd0e78/fast5_all/"
name_extension <- "gencode.v47_ENO2.tx.uf.polyA.human.noMT_protEnc.ENO2.RN7SK.ninetails_pF"

# human
groupname_hum <- "gencode.v47_ENO2.tx.uf.polyA.human.noMT_protEnc"
group_id_hum <- "human"
do_ninetails_runs_0(base_path, basename, groupname, seqsum_path, workspace_path)

# viral
groupname_vir <- NULL
group_id_vir <- NULL

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_hum,
                       groupname_2 = groupname_vir,
                       group_id_1 = group_id_hum,
                       group_id_2 = group_id_vir,
                       name_extension = name_extension)


### IVT RN7SKpolyA-002 #########################################################

ninetails_run <- groupname_hum <- groupname_vir <- group_id_hum <- group_id_vir <- NULL
base_path <- "../data/"
basename <- "RN7SKpolyA-002.guppy.hac.6.1.7"
seqsum_path <- "../input/RN7SKpolyA-002-hac.v6.1.7/"
workspace_path <- "../input/RN7SKpolyA-002-hac.v6.1.7/workspace/20250918_1051_MN32174_FAV72175_a0bce578/fast5_all"
name_extension <- ".RN7SK.tx.uf.ninetails_pF"

# ivt
groupname_ivt <- "RN7SK.tx.uf.polyA"
group_id_ivt <- "IVT_RN7SK"
do_ninetails_runs_0(base_path, basename, groupname, seqsum_path, workspace_path)

#
groupname_2 <- NULL
group_id_2 <- NULL

postproscess_ninetails(base_path = base_path,
                       basename = basename,
                       groupname_1 = groupname_ivt,
                       groupname_2 = groupname_2,
                       group_id_1 = group_id_ivt,
                       group_id_2 = group_id_2,
                       name_extension = name_extension)


###############################################################################/
