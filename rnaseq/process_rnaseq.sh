#!/bin/bash

#this is less of a script and more of a collection of the commands I used to do RNA-Seq
#I did intend for this to be a script but because I didn't have room on the cluster
#I ran things locally. May or may not have been a good idea. Anyway since I ran things locally I gave my computer plenty of time
#to cool down between commands (especially alignment) hence why I didn't run this as a script.

hisat-3n-build --base-change A,G Nvec200.fasta Nvec200

cd N1
fastqc N1_1.fq.gz
fastqc N1_2.fq.gz

trimmomatic -Xmx8g PE -threads 4 N1_1.fq.gz N1_2.fq.gz N1_forward_paired.fq.gz N1_forward_unpaired.fq.gz N1_reverse_paired.fq.gz N1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
fastqc N1_forward_paired.fq.gz
fastqc N1_reverse_paired.fq.gz

echo "Aligning N1 (this may take a while...)"
#STAR --genomeDir ../star_index --runThreadN 4 --readFilesIn N1_forward_paired.fq.gz N1_reverse_paired.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix N1
hisat-3n -t -p 4 -x ../Nvec200 -1 N1_forward_paired.fq.gz -2 N1_reverse_paired.fq.gz -S N1.sam --base-change A,G
samtools view -b -S N1.sam | samtools sort -@ 4 -o N1.sorted.bam -
samtools index N1.sorted.bam

echo "Quality check for N1"
#samtools view -H N1.sorted.bam | grep '^@RG' 
samtools flagstat N1.sorted.bam

echo "MAPQ quality check"
samtools view N1.sorted.bam | cut -f5 | sort -n | uniq -c
#actual filtering here: note -q 59, which is consistent with HISAT's unique read score
#if using STAR use -q 255
samtools view -b -F 1804 -f 2 -q 59 N1.sorted.bam > N1.filt.bam
samtools sort -o N1.filt.sort.bam N1.filt.bam
samtools index N1.filt.sort.bam

cd ../N2
fastqc N2_1.fq.gz
fastqc N2_2.fq.gz

trimmomatic -Xmx8g PE -threads 4 N2_1.fq.gz N2_2.fq.gz N2_forward_paired.fq.gz N2_forward_unpaired.fq.gz N2_reverse_paired.fq.gz N2_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
fastqc N2_forward_paired.fq.gz
fastqc N2_reverse_paired.fq.gz

echo "Aligning N2 (this may take a while...)"
hisat-3n -t -p 4 -x ../Nvec200 -1 N2_forward_paired.fq.gz -2 N2_reverse_paired.fq.gz -S N2.sam --base-change A,G
#STAR --genomeDir ../star_index --runThreadN 4 --readFilesIn N2_forward_paired.fq.gz N2_reverse_paired.fq.gz  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix N2
samtools view -b -S N2.sam | samtools sort -@ 4 -o N2.sorted.bam -
samtools index N2.sorted.bam

echo "Quality check for N2"
#samtools view -H N2.sorted.bam | grep '^@RG' 
samtools flagstat N2.sorted.bam

echo "MAPQ quality check"
samtools view N2.sorted.bam | cut -f5 | sort -n | uniq -c
samtools view -b -F 1804 -f 2 -q 59 N2.sorted.bam > N2.filt.bam
samtools sort -o N2.filt.sort.bam N2.filt.bam
samtools index N2.filt.sort.bam

cd ../N3
fastqc N3_1.fq.gz
fastqc N3_2.fq.gz

trimmomatic -Xmx8g PE -threads 4 N3_1.fq.gz N3_2.fq.gz N3_forward_paired.fq.gz N3_forward_unpaired.fq.gz N3_reverse_paired.fq.gz N3_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
fastqc N3_forward_paired.fq.gz
fastqc N3_reverse_paired.fq.gz

echo "Aligning N3 (this may take a while...)"
hisat-3n -t -p 4 -x ../Nvec200 -1 N3_forward_paired.fq.gz -2 N3_reverse_paired.fq.gz -S N3.sam --base-change A,G
#STAR --genomeDir ../star_index --runThreadN 4 --readFilesIn N3_forward_paired.fq.gz N3_reverse_paired.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix N3
samtools view -b -S N3.sam | samtools sort -@ 4 -o N3.sorted.bam -
samtools index N3.sorted.bam

echo "Quality check for N3"
#samtools view -H N3.sorted.bam | grep '^@RG' 
samtools flagstat N3.sorted.bam

echo "MAPQ quality check"
samtools view N3.sorted.bam | cut -f5 | sort -n | uniq -c
samtools view -b -F 1804 -f 2 -q 59 N3.sorted.bam > N3.filt.bam
samtools sort -o N3.filt.sort.bam N3.filt.bam
samtools index N3.filt.sort.bam

cd ..
featureCounts -p --countReadPairs -s 2 -a NV2g.20240221.gtf -o rna_counts.txt N1/N1.filt.sort.bam N2/N2.filt.sort.bam N3/N3.filt.sort.bam