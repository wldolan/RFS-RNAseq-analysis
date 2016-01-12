 #!/bin/bash

#FILENAME: bamtosam.sh

#qsub -q biochem -l nodes=1:ppn=20 -l walltime=10:00:00 bamtosam.sh

module load bioinfo
module load samtools/1.2



for file in /scratch/snyder/w/wsoltau/bamfiles/*.bam

do



##sort bam file
##samtools sort [-l level] [-m maxMem] [-o out.bam] [-O format] [-n] -T out.prefix [-@ threads] [in.bam]
samtools sort -o ${file%.*}sorted.bam -T ${file%.*}sortemp.prefix ${file}

##index bam file
###samtools index [-bc] [-m INT] aln.bam|aln.cram
samtools index ${file%.*}sorted.bam


##convert bam to sam
##samtools view [options] in.bam|in.sam|in.cram [region...]
samtools view -o ${file%.*}.sam ${file%.*}sorted.bam


done
