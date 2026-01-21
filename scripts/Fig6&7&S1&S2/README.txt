This folder contains the scripts to generate figure 6/7/S1/S2, as well as the required files to do so (mainly direclty necessary files are provided, most intermediate files in the pipeline are not provided.) (Large files need to be downloaded from ENA/zenodo first, for links see the publication.)

Steps to just create the plots:
- Install Ninetails v1.0.3 and other required packages (see publication) and run ./figures_6_7.R (required intermediate input files need to be provided in ./data/ (download from zenodo)).


Steps to run thie complete pipeline (basecalling, aligning, poly(A)-tail length estimation via nanopolish, applying ninetails, generating plots (fig 6/7/S1/S2)):

- (Note: in the HPC environment, modules as specified in the publication were used)
- Install Ninetails v1.0.3 and other required packages/ create therequired conda environment (see publication and Ninetails' installation guide (GitHub Wiki))
- Download the raw FAST5 files from zenodo and put them in their respective folders at ./raw/
- Run basecalling_and_align.sh
- (Note: FASTA files which include human references are not included in this repository for space reasons and would need to be created first)
	- gencode.v47_TB40-BAC.v1.3_ENO2.transcripts is merged from the TB40 FASTA and ENO2 FASTA (both in this repository) and the human gencode v47 transcriptome reference (to be downloaded)
	- gencode.v47_ENO2.tx is merged from the ENO2 FASTA and the human gencode v47 transcriptome reference
- Put the guppy-outputted FAST5 files in one shared folder for each dataset (./input/*/workspace/*/fast5_*/)
- Run nanopolish.sh two seperate times (once with section A for TB40 24h/48h/72h-1, normally; and once with section B for all other datasets via a slurm array (with values 1-2,5-9) (this (more complicated approach) is technically not necessary, but recomended as this way the outputs match with the rest of the pipeline. Running this only once with identical BAM file treatment is possible but would require manual changing of some scripts downstream).
- Split the nanopolish output file into human/viral/(eno2/rn7sk)-reads-only containing subfiles (adapt the names in the following awk command where necessary (to change: 'input' and 'name1'/'name2'/'name3')):
	- awk 'BEGIN {name1 = ".human.tsv"; name2 = ".hcmv.tsv"; name3 = ".eno2.tsv"; file = gensub(/\.tsv$/, "", "g", ARGV[1])} FNR==1 {print $0 > file name1; print $0 > file name2; print $0 > file name3; next} {if($2 ~ /^ENST/){print $0 >> file name1; next}; if($2 ~ /^[m]?RNA/){print $0 >> file name2; next}; if($2 == "ENO2"){print $0 >> file name3}}' input.polyA.tsv
- For the human subset: Exclude mitochondrial RNAs and only include protein coding RNAs (set the names respectively (to change: 'XXX')):
	- cat XXX.tsv | awk '{if(NR==1){print; next}; if ($2 !~ /\|MT-/ && $2 ~ /(protein_coding\|)$/){print}}' - > XXX.noMT_protEnc.tsv
	- in case of the Mock sample, the RN7SK spike-in needs to be removed (adapt the names where necessary ('name1'/'name2'/'name3')):
		- awk 'BEGIN {name1 = ".RN7SK.tsv"; name2 = ".human.tsv"; name3 = ".ENO2.tsv"; file = gensub(/\.tsv$/, "", "g", ARGV[1])} FNR==1 {print $0 > file name1; print $0 > file name2; print $0 > file name3; next} {if($2 ~ /.+RN7SK-201.+(snRNA\|)$/){print $0 >> file name1; next}; if($2 ~ /^ENST/){print $0 >> file name2; next}; if($2 == "ENO2"){print $0 >> file name3}}' NHDFpolyA-RN7SK-002.guppy.hac.6.1.7.gencode.v47_ENO2.tx.polyA.tsv
- R section (adjust names in scripts (e.g. specific environment names) where necessary):
	- Install R packages based on specification noted in the publication (if necessary, follow the installation guide on the Ninetails GitHub page)
	- Create a conda env based on the r-ninetails.yaml provided here (./r-ninetails.yaml)
	- Run fix_multiple_basecallgroups_ninetails.R [this fix in no longer necessary with Ninetails >=v1.0.56]
	- Run running_ninetails.R
	- Run figures_6_7.R
	- Plots are saved at ./plots
