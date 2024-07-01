#!/bin/bash

#SBATCH --job-name=mergeBam                # replace name
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=1024                      # replace with value for your job
#SBATCH --time=2:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add Anaconda3/2022.05
source activate HISAT2
module add SAMtools/1.16.1-GCC-11.3.0
module add BEDTools/2.30.0-GCC-11.3.0

GENOME2="~/annotations/Calbicans_A22/C_albicans_SC5314_A22" # version A22

cd ~/mapping/bam/

#for safety reasons
shopt -s nullglob
#create a bash array with all files
declare -a arr
arr=("30C_" "39C_" "PRP19_30C_" "PRP19_39C_" "SN95_30C_" "SN95_39C_" "S_")

COND=${arr[$SLURM_ARRAY_TASK_ID]}
FILENAME1="replicates/"$COND"1_A22.bam"
FILENAME2="replicates/"$COND"2_A22.bam"
FILENAME3="replicates/"$COND"3_A22.bam"

echo "$COND"

# merge bam files
samtools merge -o "$COND"pool_A22.bam $FILENAME1 $FILENAME2 $FILENAME3

# extract strand-specific reads
samtools view -b -f 0x10 "$COND"pool_A22.bam > "$COND"pool_A22_r.bam # need to change anything with paired-end data?
samtools view -b -F 0x10 "$COND"pool_A22.bam > "$COND"pool_A22_f.bam
                                
samtools depth "$COND"pool_A22_r.bam > "$COND"pool_A22_r.wig
samtools depth "$COND"pool_A22_f.bam > "$COND"pool_A22_f.wig
                                                        
cat "$COND"pool_A22_r.wig | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > "$COND"pool_A22_r.bedgraph
cat "$COND"pool_A22_f.wig | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > "$COND"pool_A22_f.bedgraph
                                                                        
rm "$COND"pool_A22_[fr].wig

igvtools count -w 1 --strands 'first' "$COND"pool_A22.bam "$COND"pool_A22.tdf "$GENOME2".genome    
		
# cleanup / compress files
rm "$COND"pool_A22_[fr].bam

gzip "$COND"pool_A22_[fr].bedgraph