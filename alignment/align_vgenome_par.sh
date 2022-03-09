#!/bin/bash

#Submits alignment for each barcode file in parallel--PREFERRED

module load fastx_toolkit
module load bowtie
module load samtools
module load bedtools2

index2='/path/to/indices' # directory to genome index 
ifqfile='/path/to/fastq/file(s)' # directory containing the fastq file(s) to map
baseoutput='/output/directory' 


cd $ifqfile
find . -maxdepth 1 -name "*BC*" -print | while read file
do
i=$(sed -e 's#.*BC\(\)#\1#' <<< "$file" | cut -f 1 -d '.')
fqfile=$ifqfile'BC'$i'.fq'
basepath=$(basename $fqfile .fq)'-'$(basename $index2)
output=$baseoutput'/'$(basename $index2)'_ALIGNMENT-Riboprof/'$(basename $fqfile .fq)'/'
echo $output

echo "module load bowtie; if [ ! -d $output ]
then
mkdir -p $output
cd $output
bowtie $index2 -q $fqfile -v 1 -m 10 -S $output$basepath'.sam';

module load samtools

samtools view -bS $output$basepath'.sam' > $output$basepath'.bam'
samtools sort $output$basepath'.bam' -o $output$basepath'_sorted.bam'
samtools index $output$basepath'_sorted.bam' $output$basepath'_sorted.bai'
samtools mpileup -BQ0 -d 500000 -f $index2'.fasta' $output$basepath'_sorted.bam' > $output$basepath'_mpileup';
python2.7 /path/to/pileup_To_counts.py $output$basepath'_mpileup' > $output$basepath'_counts';
/path/to/repair_counts.pl $output$basepath'_counts' > $output$basepath'_counts_repaired';
fi" 

done