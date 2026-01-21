This folder is designated to contain the files necessary to run ../running_ninetails.R:
- the FAST5 files (outputted by guppy (containing the basecall group), not the original sequencer output) and
- the guppy sequencing summaries
To just create the visualisations, only figures_6_7.R is necessary.

The folder structure is already set and can be downloaded from zenodo together with the guppy sequencing summary files. The FAST5 files can be created from the pipeline and put into the folder of their respective dataset (workspace/*/fast5_*/).

(If you also want to retrace the basecalling, alignment, filtering, and nanopolish step, put the original raw FAST5 files (ENA) in ./raw/ and use basecalling_and_align.sh and nanopolish.sh. as described in ../README.txt.)

