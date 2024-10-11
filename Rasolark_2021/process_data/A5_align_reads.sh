#!/bin/bash -l

#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J align_reads
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Input parameters
IND=$1

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT/A5_align_reads;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;

### Load modules
module load bioinfo-tools bwa/0.7.17 \
samtools/1.14 picard/2.23.4;

### Set up file info
find $WORKDIR/A2_quality_trim -name "*$IND*_paired.fq.gz" > $OUTDIR/samples/$IND\_1.txt;
cat $OUTDIR/samples/$IND\_1.txt | awk -F"_" '{NF-=2; gsub(" ", "_"); print}' | \
sort -u > $OUTDIR/samples/$IND\_2.txt;
cat $OUTDIR/samples/$IND\_2.txt | rev | cut -d'/' -f1 | rev \
> $OUTDIR/samples/$IND\_3.txt;
paste -d '\t' $OUTDIR/samples/$IND\_2.txt \
$OUTDIR/samples/$IND\_3.txt > $OUTDIR/samples/$IND.txt;

### Align reads to reference and add @RG header
for FILES in $(seq 1 1 $(cat $OUTDIR/samples/$IND.txt | wc -l)); do

# Input parameters
FILE=$(cat $OUTDIR/samples/$IND.txt | head -n $FILES | tail -n 1 | cut -f1);
ID=$(cat $OUTDIR/samples/$IND.txt | head -n $FILES| tail -n 1 | cut -f2);
LIB=$(cat $METADATA/sample_info.txt | awk -F'\t' -v IND=$IND '$1==IND {print $6}');
LANE=$(zcat $(find $WORKDIR/A2_quality_trim -name "$ID*R1_paired.fq.gz") | \
head -n1 | cut -d' ' -f1 | cut -d'@' -f2 | cut -d':' -f1,2,3,4);
PU=$(zcat $(find $WORKDIR/A2_quality_trim -name "$ID*R1_paired.fq.gz") \
| head -n1 | cut -d' ' -f1 | cut -d'@' -f2 | cut -d':' -f3);

# Align paired reads
bwa mem -t 20 -M $REF $FILE\_R1_paired.fq.gz $FILE\_R2_paired.fq.gz \
-R "@RG\tID:$IND-$LANE\tSM:$IND\tLB:$LIB\tPU:$PU\tPL:illumina" | \
samtools view -b -@ 20 > $OUTDIR/temp1/$ID.bam;

# Sort bam file
samtools sort -n -T $OUTDIR/temp1/$ID -@ 20 \
$OUTDIR/temp1/$ID.bam > $OUTDIR/temp2/$ID\_qsorted.bam;
done;

# Merge bams from the same individual if more than one bam file
if [ $(cat $OUTDIR/samples/$IND.txt | wc -l) -gt 1 ]; then
  BAMS=$(find $OUTDIR/temp2 -name "*-$IND\_*");
  samtools merge $OUTDIR/temp3/$IND\_qsorted_merged.bam -n -@ 20 $BAMS;
else
  mv $OUTDIR/temp2/$ID\_qsorted.bam $OUTDIR/temp3/$IND\_qsorted_merged.bam;
fi;

### Remove duplicates, and index bam files

# Input parameters
PCR_FREE=$(cat $METADATA/sample_info.txt | awk -F'\t' -v IND=$IND '$1==IND {print $5}');
if [ $PCR_FREE == "Yes" ]; then
  DUP=$(echo REMOVE_SEQUENCING_DUPLICATES=true);
else
  DUP=$(echo REMOVE_DUPLICATES=true);
fi;

# Mark and remove duplicates
java -Xmx100g -jar $PICARD_ROOT/picard.jar MarkDuplicates \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
$DUP ASSUME_SORT_ORDER=queryname \
I=$OUTDIR/temp3/$IND\_qsorted_merged.bam \
O=$OUTDIR/temp4/$IND\_qsorted_merged_nodup.bam \
M=$OUTDIR/duplicates_data/$IND\_duplicate_data.txt \
R=$REF TMP_DIR=$OUTDIR/temp5;

# Filter and sort bam file by coordinate
samtools view -f2 -F260 -q20 -b -@ 20 \
$OUTDIR/temp4/$IND\_qsorted_merged_nodup.bam | \
samtools sort -T $OUTDIR/$IND.bam -@ 20 \
> $OUTDIR/$IND\_qsorted_merged_nodup_filtered_sorted.bam;
samtools index $OUTDIR/$IND\_qsorted_merged_nodup_filtered_sorted.bam -@ 20;
