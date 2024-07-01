#!/bin/bash

#SBATCH --job-name=pipeline                # replace name
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=4096                      # replace with value for your job
#SBATCH --time=10:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job


pwd

#for safety reasons
shopt -s nullglob
#create a bash array with all files
VAR=$(ls 2[23]*/*_R1_*.fastq.gz | grep -v "PMGC_Control_Cowen_S28_L001_R1_001.fastq.gz")

declare -a arr
arr=($VAR)

R1=${arr[$SLURM_ARRAY_TASK_ID]}
R2=$(echo $R1 | sed 's/_R1_/_R2_/')
OUTNAME=$(echo $R1 | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/_S.*_L001_R1_.*.fastq.gz//')

echo $R2
echo $OUTNAME

sh ~/scripts/adaptor_trim_20230623.sh $R1 $R2 $OUTNAME

mv *_R[12]t.fastq.gz fastq_trim/
mv QC_$OUTNAME fastq_trim/

cd mapping
sh ~/scripts/mapping_HISAT2_Ca_20230623.sh ../fastq_trim/"$OUTNAME"_R1t.fastq.gz ../fastq_trim/"$OUTNAME"_R2t.fastq.gz $OUTNAME
