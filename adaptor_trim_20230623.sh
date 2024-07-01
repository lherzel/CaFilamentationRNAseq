#!/bin/bash

module add cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2
module add FastQC/0.11.9-Java-11

R1=$1
R2=$2
OUTNAME=$3

cutadapt --cores=0 -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT -e 0.1 -O 10 -m 16 --output="$OUTNAME"_R1t.fastq.gz --paired-output="$OUTNAME"_R2t.fastq.gz $R1 $R2

mkdir QC_$OUTNAME
fastqc -o QC_$OUTNAME $R1 $R2 "$OUTNAME"_R1t.fastq.gz "$OUTNAME"_R2t.fastq.gz

