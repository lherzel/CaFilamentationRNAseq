#! /bin/sh

BAM=$1
ANNO=$2
OUT=$3
READLENGTH="101"

#Use intron annotation
#Overlap reads with introns
#Classify reads
#	Junction unsplit
#	Junction split (can be any, even completely spanning)
#	Internal i
#Normalize internal counts by length L in nt
#	i / L, e.g. 100/1000 - 0.1
#IR score
#(Jct unsplit + norm internal) / (all jct + norm internal)
#FOLDERNAME="$OUT"_IRdetector_"$RANDOM"
FOLDERNAME="$OUT"_CaIRdetector

mkdir $FOLDERNAME
cd $FOLDERNAME

# split bam-file into individual reads for pairs
#samtools view -F 4 $BAM | wc -l > counts # read counts

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

cat ol1 ol2 | awk ' $NF>=3 { print $0 } ' > ol
cat ol | awk ' $10==1 { print $0 }' > unspliced
cat unspliced | awk ' { if ($2>=$14 && $3<=$15) print $0"\tintronic"; else print $0"\tunspliced" } ' > unsplicedAdj

cat ol | awk ' $10>1 { print $0 }' > spliced
cat spliced | awk ' { print $11 }' | sed 's/,/\t/g' > blocks.read

paste spliced blocks.read > jct

cat jct | awk ' $10==2 { print $0 } ' > jct1
#cat jct | awk ' $10==3 { print $0 } ' > jct2
#cat jct | awk ' $10==4 { print $0 } ' > jct3

# if [ -s jct3 ]
# 	then
# 		echo "junctions with 3 introns exist"
# 	else
# 		rm jct3
# 	fi 
	
cat jct1 | awk ' { 
	if (($2+$20)==$14 && ($3-$21)==$15) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\tspliced";
	else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\talternative" 
	} ' > jct1Adj

# if read spans multiple introns check that at least 3 nt in all blocks
# do not include multi-intron spanning reads for now... too complicated
# cat jct2 | awk ' $20>=3 && $22>=3 { print $0 } ' > jct2Adj
# cat jct2 | awk ' { if ($20<3) print $21"\t"$22"\t"$15"\t"$3-$22"\t"$15; if ($22<3) print $20"\t"$21"\t"$2+$20"\t"$14; else print $20"\t"$21"\t"$22 } ' | head

cat unsplicedAdj jct1Adj | awk ' { print $1"\t"$2"\t"$3"\t"$4"\t"$10"\t"$11"\t"$12"\t"$14"\t"$15"\t"$16"\t"$18"\t"$19"\t"$20 } ' > classified.reads

Rscript ~/scripts/CaIR_detector_LHerzel_20230809.R $READLENGTH
# outputs quantified file

mv quantified "$OUT"_IRcounts.txt

#rm classified.reads counts1 counts2
#rm r1.bam r2.bam ol blocks.read spliced unspliced unsplicedAdj jct1Adj ol1 ol2
