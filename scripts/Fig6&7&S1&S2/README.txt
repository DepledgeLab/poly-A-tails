Steps to run this pipeline (basecall, align, poly(A)-tail estimation via nanopolish, apply ninetails, generate plots (fig 6/7/S1/S2)):

- (Note: in the HPC environment, modules as specified in the paper were used)
- Put the raw FAST5 and fastq files in their respective folders at ../input (full relative paths needed are noted in the 'dataset_info_*.txt' files; For TB40 24h/48h/72h (1) paths are noted in basecalling_and_align.sh)
- Run basecalling_and_align.sh
- (Note: fasta files which include human references are not aincluded in this repo for storage reasons)
	- gencode.v47_TB40-BAC.v1.3_ENO2.transcripts is merged from the TB40 fasta and ENO2 fasta (both in this repo) and the human gencode v47 transcriptome reference
	- gencode.v47_ENO2.tx does not include TB40
- Run nanopolish.sh (Note: For HCMV (TB40) 25h/48h/72h (1), use the commented out sections; for the other datasets use the slurm array with the 'dataset_info_*.txt' files (1-2,5-9))
- Split the nanopolish output file into human/viral/(eno2/rn7sk) (adapt names where necessary (input and name1/2/3)):
	- awk 'BEGIN {name1 = ".human.tsv"; name2 = ".hcmv.tsv"; name3 = ".eno2.tsv"; file = gensub(/\.tsv$/, "", "g", ARGV[1])} FNR==1 {print $0 > file name1; print $0 > file name2; print $0 > file name3; next} {if($2 ~ /^ENST/){print $0 >> file name1; next}; if($2 ~ /^[m]?RNA/){print $0 >> file name2; next}; if($2 == "ENO2"){print $0 >> file name3}}' input.polyA.tsv
- For human subset: Exclude mitochondrial RNAs and only include protein coding RNAs (set the names respectively ('XXX')):
	- cat XXX.tsv | awk '{if(NR==1){print; next}; if ($2 !~ /\|MT-/ && $2 ~ /(protein_coding\|)$/){print}}' - > XXX.noMT_protEnc.tsv
	- in case of the Mock sample, the RN7SK spike-in needs to be removed:
		- awk 'BEGIN {name1 = ".RN7SK.tsv"; name2 = ".human.tsv"; name3 = ".ENO2.tsv"; file = gensub(/\.tsv$/, "", "g", ARGV[1])} FNR==1 {print $0 > file name1; print $0 > file name2; print $0 > file name3; next} {if($2 ~ /.+RN7SK-201.+(snRNA\|)$/){print $0 >> file name1; next}; if($2 ~ /^ENST/){print $0 >> file name2; next}; if($2 == "ENO2"){print $0 >> file name3}}' NHDFpolyA-RN7SK-002.guppy.hac.6.1.7.gencode.v47_ENO2.tx.polyA.tsv
- R section (adjust names in scripts (e.g. specific env names) were necessary):
	- Install R packages based on specification noted in the paper (if necessary, follow the installation guide on the ninetails GitHub page)
	- Create conda env based on r-ninetails.yaml
	- Run fix_multiple_basecallgroups_ninetails.R
	- Run running_ninetails.R
	- Run figures_6_7.R

