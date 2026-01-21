#!/bin/bash
#SBATCH --job-name=basecalling_and_align  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=01:00:00 # Time limit hrs:min:sec
#SBATCH --output=../../logs/%x_%j.log # Standard output and error log
#SBATCH -p <partition>

# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> Basecalling ONT DRS datasets (hac, guppy) and Alignment (minimap2) to
# the reference transcriptomes
# (In addition to variables) needs following files/directories:
#   - Hybrid transcriptome reference (e.g. Human + Viral + ENO2) [.fasta]
#   - raw sequencing data [.fast5 direc]

# /================================================================================================\
# |                                           VARIABLES                                            |
# |________________________________________________________________________________________________|

OUT=./data
IN=./raw
GUPPY_PATH=/path/to/ont-guppy-6.1.7/bin/guppy_basecaller
GUPPY_VER=6.1.7
GUPPY_MOD=hac
BASECALL_DIREC_BASE=./input

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

module load SAMtools
module load minimap2
module load HDF5
module load VBZ-Compression

#  ------------------------------------------------------------------------------------------------

    echo -e "\n"
    echo "Modules loaded successfully."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                             CODE                                               |
# |________________________________________________________________________________________________|

while read NAME_BASE TX_REF TX_ABBR LOREM IPSUM;

do

BASECALL_DIREC="$BASECALL_DIREC_BASE"/"$NAME_BASE"-"$GUPPY_MOD".v"$GUPPY_VER"
NAME="$NAME_BASE".guppy."$GUPPY_MOD"."$GUPPY_VER"

        echo -e "\n================\n Begun "$NAME" \n================\n\n"

### Basecalling
        echo -e "\n================\n Begun basecalling \n================\n\n"
        SECONDS=0

        $GUPPY_PATH -i $IN/$NAME_BASE/ -s $BASECALL_DIREC -c rna_r9.4.1_70bps_hac.cfg -r --calib_detect --trim_strategy rna --reverse_sequence true -x auto --fast5_out

        duration=$SECONDS
        echo -e "\nFinished basecalling after $((duration / 3600)) hour(s), $(((duration / 60) % 60)) minute(s) and $((duration % 60)) second(s).\n----------\n"

### Processing
        echo -e "\n================\n Begun processing \n================\n\n"
        SECONDS=0

        cat $BASECALL_DIREC/pass/*.fastq > $OUT/"$NAME".fastq
        samtools faidx $TX_REF -o "$TX_REF".fai

        duration=$SECONDS
        echo -e "\nFinished pre-processing after $((duration / 3600)) hour(s), $(((duration / 60) % 60)) minute(s) and $((duration % 60)) second(s).\n----------\n"

### Alignment
        echo -e "\n================\n Begun aligning \n================\n\n"
        SECONDS=0

        minimap2 -t 8 -ax map-ont -L -p 0.99 -uf $TX_REF $OUT/"$NAME".fastq > $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.sam
        samtools view -F2324 -b -o $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.bam $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.sam
        samtools sort -o $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.sorted.bam $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.bam
        samtools index $OUT/revised_with_aligning_uf/"$NAME"."$TX_ABBR".uf.sorted.bam

        duration=$SECONDS
        echo -e "\nFinished aligning after $((duration / 3600)) hours, $(((duration / 60) % 60)) minutes and $((duration % 60)) seconds.\n----------\n"


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

