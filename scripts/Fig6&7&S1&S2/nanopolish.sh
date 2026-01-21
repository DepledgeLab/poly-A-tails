#!/bin/bash
#SBATCH --job-name=nanopolish  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=../../logs/%x_%j.log # Standard output and error log
#SBATCH -p <partition>

# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> Poly(A)-tail estimation with nanopolish.

# (In addition to variables) needs following files/directories:
#   - Hybrid transcriptome reference (e.g. Human + Viral + ENO2 (spike-in)) [.fasta])
#   - sequencing data (Guppy output) [.fast5 direc]
#   - Guppy basecalled data [.fastq direc]
#   - Guppy sequence summary [.txt]

# Run this script 2 seperate times. To apply nanopolish to the TB40 24h/48h/72h-1 datasets, use
# section A (and comment out section B. The useed paths are included in this script).
# To apply nanopolish to the other datasets, use section B (and comment out section A) via a slurm 
# array. The used paths for each dataset are in the array-loaded files (./dataset_info_*.txt).

# /================================================================================================\
# |                                           VARIABLES                                            |
# |________________________________________________________________________________________________|

OUT=./data
GUPPY_VER=6.1.7
GUPPY_MOD=hac
NANOPOLISH_PATH=/path/to/nanopolish


# TB40
REF_TX1=../../transcriptomes/GRCh38.p13.gencode.v42.transcripts.fasta
REF_TX2=../../transcriptomes/TB40-BAC-tx-v1.3.fasta
SET1=GRCh38.tx
SET2=TB40txome

# TB40 72h-1
BASECALL_DIREC_72=./input/TB40_72h-hac.v6.1.7
NAME_72=TB40_72h
FAST5_DIREC_72=./input/TB40_72h-hac.v6.1.7/workspace/20200207_1844_MN24978_FAL81867_7f514654/fast5_all

# TB40 48h
BASECALL_DIREC_48=./input/TB40_48h-hac.v6.1.7
NAME_48=TB40_48h
FAST5_DIREC_48=./input/TB40_48h-hac.v6.1.7/workspace/20200130_1934_MN24978_FAL86847_be4e2264/fast5_all

# TB40 24h
BASECALL_DIREC_24=./input/TB40_24h-hac.v6.1.7
NAME_24=TB40_24h
FAST5_DIREC_24=./input/TB40_24h-hac.v6.1.7/workspace/fast5_pass


# /================================================================================================\
# |                                         INITIALISATION                                         |
# |________________________________________________________________________________________________|

    starttime=$(date +"%D %T")
    echo -e "\n"
    echo -e "           ============================================"
    echo -e "           =     started at $starttime         ="
    echo -e "           ============================================"

    echo "<<< Start Job '$SLURM_JOB_NAME' (Job ID: $SLURM_JOB_ID) on $HOSTNAME by $SLURM_JOB_USER >>>"  # Display job start information
    echo "<<< Using $SLURM_CPUS_ON_NODE CPU(s). Running with $SLURM_MEM_PER_NODE MiB and $SLURM_NTASKS task(s). >>>" # Display job execution information
    echo "Script successfully started."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                            MODULES                                             |
# |________________________________________________________________________________________________|

module load HDF5
module load VBZ-Compression
module load SAMtools

#  ------------------------------------------------------------------------------------------------

    echo -e "\n"
    echo "Modules loaded successfully."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                             CODE                                               |
# |________________________________________________________________________________________________|

# Run either section A or section B, comment out the other

# Section A: TB40 24h/48h/72h-1

	echo -e "\n================\n Begun "$NAME_24"\n================\n\n"
### Nanopolish-index
	SECONDS=0

	$NANOPOLISH_PATH index --directory=$FAST5_DIREC_24/ --sequencing-summary=$BASECALL_DIREC_24/sequencing_summary.txt $OUT/"$NAME_24".fastq

        duration=$SECONDS
        echo -e "\nFinished nanopolish indexing after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"

### Nanopolish-polya
        SECONDS=0
	
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_24".fastq --bam=$OUT/$NAME_24.$SET1.sorted.bam --genome=$REF_TX1 > $OUT/$NAME_24.$SET1.polyA.tsv
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_24".fastq --bam=$OUT/$NAME_24.$SET2.sorted.bam --genome=$REF_TX2 > $OUT/$NAME_24.$SET2.polyA.tsv

        duration=$SECONDS
        echo -e "\nFinished Nanopolish-polya after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"


	echo -e "\n================\n Begun "$NAME_48"\n================\n\n"
### Nanopolish-index
	SECONDS=0

	$NANOPOLISH_PATH index --directory=$FAST5_DIREC_48/ --sequencing-summary=$BASECALL_DIREC_48/sequencing_summary.txt $OUT/"$NAME_48".fastq

        duration=$SECONDS
        echo -e "\nFinished nanopolish indexing after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"

### Nanopolish-polya
        SECONDS=0
	
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_48".fastq --bam=$OUT/$NAME_48.$SET1.sorted.bam --genome=$REF_TX1 > $OUT/$NAME_48.$SET1.polyA.tsv
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_48".fastq --bam=$OUT/$NAME_48.$SET2.sorted.bam --genome=$REF_TX2 > $OUT/$NAME_48.$SET2.polyA.tsv

        duration=$SECONDS
        echo -e "\nFinished Nanopolish-polya after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"


	echo -e "\n================\n Begun "$NAME_72"\n================\n\n"
### Nanopolish-index
	SECONDS=0

	$NANOPOLISH_PATH index --directory=$FAST5_DIREC_72/ --sequencing-summary=$BASECALL_DIREC_72/sequencing_summary.txt $OUT/"$NAME_72".fastq

        duration=$SECONDS
        echo -e "\nFinished nanopolish indexing after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"

### Nanopolish-polya
        SECONDS=0
	
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_72".fastq --bam=$OUT/$NAME_72.$SET1.sorted.bam --genome=$REF_TX1 > $OUT/$NAME_72.$SET1.polyA.tsv
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME_72".fastq --bam=$OUT/$NAME_72.$SET2.sorted.bam --genome=$REF_TX2 > $OUT/$NAME_72.$SET2.polyA.tsv

        duration=$SECONDS
        echo -e "\nFinished Nanopolish-polya after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"



# Section B: All other datasets

while read NAME_BASE REF_TX ABBR BASECALL_DIREC FAST5_DIREC;

do

NAME="$NAME_BASE".guppy."$GUPPY_MOD"."$GUPPY_VER"


	echo -e "\n================\n Begun "$NAME"\n================\n\n"
### Nanopolish-index
	SECONDS=0

	$NANOPOLISH_PATH index --directory=$FAST5_DIREC/ --sequencing-summary=$BASECALL_DIREC/sequencing_summary.txt $OUT/"$NAME".fastq

        duration=$SECONDS
        echo -e "\nFinished nanopolish indexing after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"


### Nanopolish-polya
        SECONDS=0
	
	$NANOPOLISH_PATH polya --threads=8 --reads=$OUT/"$NAME".fastq --bam=$OUT/$NAME.$ABBR.uf.sorted.bam --genome=$REF_TX > $OUT/$NAME.$ABBR.uf.polyA.tsv

        duration=$SECONDS
        echo -e "\nFinished Nanopolish-polya after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"

done < dataset_info_${SLURM_ARRAY_TASK_ID}.txt


# /================================================================================================\
# |                                             ENDING                                             |
# |________________________________________________________________________________________________|

    echo -e "\n"

    endtime=$(date +"%D %T")
    echo -e "\n"
    echo -e "           ============================================"
    echo -e "           =     finished at $endtime        ="
    echo -e "           ============================================"


# =================================================================================================

