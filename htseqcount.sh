#!/bin/bash

#FILENAME: htseqcount.sh

#qsub -q biochem -l walltime=10:00:00 htseqcount.sh
#qsub -q biochem -l nodes=2:ppn=20 htseqcount.sh

module load bioinfo
module load htseq/0.6.1


for file in /scratch/snyder/w/wsoltau/samfiles/*.sam

do

##htseq-count [options] <alignment_file> <gff_file>

htseq-count -s reverse -r pos -m union -t exon -i gene_id ${file} /scratch/snyder/w/wsoltau/corrected.gtf > ${file%.*}_count.txt

done
