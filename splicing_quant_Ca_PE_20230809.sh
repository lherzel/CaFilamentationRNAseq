#!/bin/bash

#SBATCH --job-name=CaSpl                # replace name
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=1024                      # replace with value for your job
#SBATCH --time=1:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add SAMtools/1.9-foss-2018b
module add BEDTools/2.27.1-foss-2018b
module add R/3.5.1-foss-2018b

LC_CTYPE="en_US.UTF-8"

cd ~/IRquant/
#pwd

#for safety reasons
shopt -s nullglob
#create a bash array with all files
VAR=$(ls ~/mapping/bam/replicates/*.bam)

declare -a arr
arr=($VAR)

ANNO="~/annotations/Calbicans_A22/Ca22_spliceosome_introns_all_woComplicated.bed"
BAM=${arr[$SLURM_ARRAY_TASK_ID]}
NAME=$(echo $BAM | awk ' BEGIN { FS="/" }; { print $NF } ' | sed 's/_A22.bam//')

#echo $BAM
echo $NAME

# run IR and AS quantification script
# intersected with an R script that also in ~/scripts/CaIR_detector_LHerzel_20230809.R
~/scripts/CaIR_detector_LHerzel_20230809.sh $BAM $ANNO $NAME
