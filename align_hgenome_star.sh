#!/bin/bash

#Aligns barcode-separated reads to human genome using star--PREFERRED 

module load star

index1='/path/to/Ribosomal_indices' # directory to ribosome index 
index2='/path/to/indices' # directory to genome index 
ifqfile='/path/to/fastq/file(s)' # directory containing the fastq file(s) to map
baseoutput='/output/directory' 
gtf='/path/to/annotation.gtf'

cd $ifqfile
find . -maxdepth 1 -name "*BC*" -print | while read file
do
i=$(sed -e 's#.*BC\(\)#\1#' <<< "$file" | cut -f 1 -d '.')
fqfile=$ifqfile'BC'$i'.fq'
basepath=$(basename $fqfile .fq)'-'$(basename $index2)
echo $basepath
output=$baseoutput'/'$(basename $index2)'_ALIGNMENT-Riboprof/'$(basename $fqfile .fq)'/'

echo "module load star;if [ ! -d $output ]
then
mkdir -p $output
cd $output

STAR --genomeDir $index1 \
     --readFilesIn $fqfile \
     --outFileNamePrefix $output \
     --outFilterMismatchNoverLmax 0.04 \
     --outSAMtype None \
     --outReadsUnmapped Fastx

STAR --genomeDir $index2 \
     --readFilesIn $output'Unmapped.out.mate1' \
     --outFileNamePrefix $output$basepath \
     --outFilterMismatchNoverLmax 0.04 \
     --outSAMstrandField intronMotif \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM


featureCounts -t exon -g gene_id -a $gtf -o $output$basepath'_fcCounts.count' $output$basepath'Aligned.sortedByCoord.out.bam'
python /path/to/reduce_fc.py $output$basepath'_fcCounts.count'
"
done
