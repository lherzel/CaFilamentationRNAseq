#! /bin/sh

BAM=$1
ANNO=$2
OUT=$3

FOLDERNAME="$OUT"_geneCounts

mkdir $FOLDERNAME
cd $FOLDERNAME

# consider to include condition for single-end and paired-end data [20230830]

# first read
samtools view -h -b -f 0x40 $BAM > r1.bam
# second read
samtools view -h -b -f 0x80 $BAM > r2.bam 

# overlap with annotated introns
# for first read require different strand
# for second read require same strand
bedtools intersect -bed -a r1.bam -b $ANNO -S -wo > ol1
bedtools intersect -bed -a r2.bam -b $ANNO -s -wo > ol2

samtools view r1.bam | wc -l > counts1
samtools view r2.bam | wc -l > counts2

#cat ol1 ol2 | awk ' $NF>=3 { print $16 } ' | sort -n | uniq -c | awk ' { print $2"\t"$1 } ' > ol
awk '$NF >= 3 { count[$16]++ } END { for (key in count) print key, count[key] }' ol1 ol2 > ol

#count rows per gene ID to get counts per gene
# add total counts per genes
wc -l ol1 | awk ' { print $1 } ' > countsG1
wc -l ol2 | awk ' { print $1 } ' > countsG2
paste countsG1 countsG2 | awk ' { print "total_gene\t"$1+$2 } ' > totalG

# add total mapped reads, including intergenic
paste counts1 counts2 | awk ' { print "total_mapped\t"$1+$2 } ' > total

# add read counts to top rows
cat total totalG ol > "$OUT"_readsPerGene.txt

rm counts1 counts2 countsG1 countsG2 total totalG 
rm r1.bam r2.bam ol ol1 ol2
