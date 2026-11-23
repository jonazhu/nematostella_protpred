#!/bin/bash

#this is less of a script and more of a collection of the commands I used to do RNA-Seq
#I did intend for this to be a script but because I didn't have room on the cluster
#I ran things locally. May or may not have been a good idea. Anyway since I ran things locally I gave my computer plenty of time
#to cool down between commands (especially alignment) hence why I didn't run this as a script.

hisat-3n-build --base-change A,G jaNemVect1.1.fna NemVect_target
    
for s in N1 N2 N3; do
    echo "-----------------------------------------"
    echo "Processing ${s}"
    fastqc ${s}_1.fq.gz
    fastqc ${s}_2.fq.gz

    trimmomatic -Xmx8g PE -threads 4 ${s}_1.fq.gz ${s}_2.fq.gz ${s}_forward_paired.fq.gz ${s}_forward_unpaired.fq.gz ${s}_reverse_paired.fq.gz ${s}_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:60 HEADCROP:10
    fastqc ${s}_forward_paired.fq.gz
    fastqc ${s}_reverse_paired.fq.gz

    echo "Aligning ${s} (this may take a while...)"
    #STAR --genomeDir ../star_index --runThreadN 4 --readFilesIn N1_forward_paired.fq.gz N1_reverse_paired.fq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix N1
    hisat-3n -t -p 4 -x NemVect_target -1 ${s}_forward_paired.fq.gz -2 ${s}_reverse_paired.fq.gz --base-change A,G | samtools sort -@ 4 -o ${s}.sorted.bam -
    samtools index ${s}.sorted.bam

    echo "Quality check for ${s}"
    #samtools view -H N1.sorted.bam | grep '^@RG' 
    samtools flagstat N1.sorted.bam

    echo "MAPQ quality check"
    samtools view ${s}.sorted.bam | cut -f5 | sort -n | uniq -c
    #actual filtering here: note -q 59, which is consistent with HISAT's unique read score
    #if using STAR use -q 255
    samtools view -b -F 1804 -f 2 -q 59 ${s}.sorted.bam > ${s}.filt.bam
    samtools sort -o ${s}.filt.sort.bam ${s}.filt.bam
    samtools index ${s}.filt.sort.bam

done

featureCounts -p --countReadPairs -s 2 -a genomic.gtf -o rna_counts.txt N1.filt.sort.bam N2.filt.sort.bam N3.filt.sort.bam