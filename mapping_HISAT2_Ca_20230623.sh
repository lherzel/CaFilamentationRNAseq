#!/bin/bash

module add Anaconda3/2022.05
source activate HISAT2

# conda create --name HISAT2
# conda install -c bioconda igvtools
# conda install -c biobuilds fastx-toolkit # not needed here (yet)

module add Python/3.10.8-GCCcore-12.2.0
module add HISAT2/2.2.1-gompi-2022a
module add SAMtools/1.16.1-GCC-11.3.0
module add BEDTools/2.30.0-GCC-11.3.0

echo Today date is:
date
echo Your current directory is:
pwd

GENOME2="~/annotations/Calbicans_A22/C_albicans_SC5314_A22" # version A22

R1=$1
R2=$2
OUTNAME=$3 #$(echo $FILE | awk 'BEGIN { FS = "/" } ; { print $NF } ' | sed 's/.gz//' | sed 's/.fastq//')
echo $OUTNAME
						
hisat2 --summary-file "$OUTNAME"_A22_summary.txt \
	--max-intronlen 2000 \
	-x $GENOME2 \
	--rna-strandness FR \
	-1 $R1 -2 $R2  \
	-S "$OUTNAME"_A22.sam \
	--novel-splicesite-outfile "$OUTNAME"_novelSS.txt

# convert sam to bam and remove any unmapped reads to reduce file size		
# repeat bam, bedgraph & tdf file generation for genome version A22
samtools view -S -b "$OUTNAME"_A22.sam > "$OUTNAME"_tmp.bam
samtools view -b -F 4 "$OUTNAME"_tmp.bam > "$OUTNAME"_A22.bam
samtools sort "$OUTNAME"_A22.bam > "$OUTNAME"_tmp.bam
mv "$OUTNAME"_tmp.bam "$OUTNAME"_A22.bam
samtools index "$OUTNAME"_A22.bam      

# extract strand-specific reads
samtools view -b -f 0x10 "$OUTNAME"_A22.bam > "$OUTNAME"_A22_r.bam # need to change anything with paired-end data?
samtools view -b -F 0x10 "$OUTNAME"_A22.bam > "$OUTNAME"_A22_f.bam
                                
samtools depth "$OUTNAME"_A22_r.bam > "$OUTNAME"_A22_r.wig
samtools depth "$OUTNAME"_A22_f.bam > "$OUTNAME"_A22_f.wig
                                                        
cat "$OUTNAME"_A22_r.wig | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > "$OUTNAME"_A22_r.bedgraph
cat "$OUTNAME"_A22_f.wig | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > "$OUTNAME"_A22_f.bedgraph
                                                                        
rm "$OUTNAME"_A22_[fr].wig

igvtools count -w 1 --strands 'first' "$OUTNAME"_A22.bam "$OUTNAME"_A22.tdf "$GENOME2".genome    
		
# cleanup / compress files
rm "$OUTNAME"_A22.sam "$OUTNAME"_A22_[fr].bam

if [ $(echo $FILE | sed 's/.*.fast//') = "q" ]
	then
		gzip $FILE
	fi

gzip "$OUTNAME"_A22_[fr].bedgraph

