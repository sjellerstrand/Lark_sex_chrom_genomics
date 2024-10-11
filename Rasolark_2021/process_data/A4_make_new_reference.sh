#!/bin/bash -l

#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -J make_new_reference
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT1=Rasolark_2021;
PROJECT2=Skylark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT1;
REF=$WORKDIR2/B1_prepare_reference\
/GCA_902810485.1_skylark_genome_genomic_no_IUPAC.fasta;
METADATA=$WORKDIR1/metadata;
IND=TT95885;

### Load modules
module load bioinfo-tools bwa/0.7.17 samtools/1.14 \
freebayes/1.3.2 vcftools/0.1.16 vcflib/1.0.1 vt/0.5772 \
bcftools/1.14 python/3.9.5 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR0/A4_make_new_reference;
OUTDIR0=$OUTDIR0/A4_make_new_reference;
mkdir $OUTDIR0/bamfiles $OUTDIR0/vcf_whole_genome;

### Set up file info
OUTDIR=$OUTDIR0/bamfiles;
find $WORKDIR1/A2_quality_trim -name "*_paired.fq.gz" | grep "$IND" > $OUTDIR/INDS1.txt;
cat $OUTDIR/INDS1.txt | awk -F"_" '{NF-=2; gsub(" ", "_"); print}' | \
sort -u > $OUTDIR/INDS2.txt;
cat $OUTDIR/INDS2.txt | rev | cut -d'/' -f1 | rev \
> $OUTDIR/INDS3.txt;
paste -d '\t' $OUTDIR/INDS2.txt \
$OUTDIR/INDS3.txt > $OUTDIR/INDS.txt;

### Align reads to reference and add @RG header
export REF WORKDIR1 OUTDIR METADATA;

## Define function
align_reads() {

# Input parameters
FILE=$1;
ID=$2;
IND=$(echo $ID | cut -d'_' -f1 | rev | awk -F'-' -v OFS='-' '{NF-=2; print}' | rev);
LIB=$(cat $METADATA/sample_info.txt | awk -F'\t' -v IND=$IND '$1==IND {print $6}');
LANE=$(zcat $(find $WORKDIR1/A2_quality_trim -name "$ID*R1_paired.fq.gz") | \
head -n1 | cut -d' ' -f1 | cut -d'@' -f2 | cut -d':' -f1,2,3,4);
PU=$(zcat $(find $WORKDIR1/A2_quality_trim -name "$ID*R1_paired.fq.gz") \
| head -n1 | cut -d' ' -f1 | cut -d'@' -f2 | cut -d':' -f3);

# Align paired reads
bwa mem -t 20 -M $REF $FILE\_R1_paired.fq.gz $FILE\_R2_paired.fq.gz \
-R "@RG\tID:$IND-$LANE\tSM:$IND\tLB:$LIB\tPU:$PU\tPL:illumina" | \
samtools view -b -@ 20 > $OUTDIR/$ID.bam;

# Filter and sort bam file by coordinate
samtools view -f2 -F260 -q20 -b -@ 20 $OUTDIR/$ID.bam | \
samtools sort -T $OUTDIR/$ID.bam -@ 20 \
> $OUTDIR/$ID\_sorted.bam;
#rm $OUTDIR/$ID.bam;
};

## Excecute function in parallell
export -f align_reads;
parallel --colsep '\t' 'align_reads {1} {2}' :::: $OUTDIR/INDS.txt;

## Merge bams from the same individual and index bam files

# Input parameters
BAMS=$(find $OUTDIR/ -name "*_sorted*");

# Merge bams from the same individual
samtools merge $OUTDIR/$IND\_sorted_merged.bam -@ 20 $BAMS;
samtools index $OUTDIR/$IND\_sorted_merged.bam;

### Setup file info
OUTDIR=$OUTDIR0/vcf_whole_genome;
fasta_generate_regions.py $REF.fai 100000 > $OUTDIR/vcf_regions.txt;
regions=$OUTDIR/vcf_regions.txt;

### Call variants in parallel with freebayes
freebayes-parallel $regions 20 -f $REF \
-b $OUTDIR0/bamfiles/$IND\_sorted_merged.bam --standard-filters \
> $OUTDIR/$IND\_whole_genome.vcf;
bgzip $OUTDIR/$IND\_whole_genome.vcf;
tabix -p vcf $OUTDIR/$IND\_whole_genome.vcf.gz;

### Filter vcf and remove indels

## Define functions
FUNCTIONS=$MAINDIR/scripts/$PROJECT1;
correct_halfcalls=$FUNCTIONS/process_data/correct_halfcalls.py;

## Filter 1
vcftools --gzvcf $OUTDIR/$IND\_whole_genome.vcf.gz \
--exclude-bed $WORKDIR2/B1_prepare_reference/$PROJECT2\_masked_repeats.bed \
--minDP 5 \
--minQ 30 \
--recode --recode-INFO-all --stdout | \
vcffilter \
-f "QA = 30 | QA > 30" \
-f "MQMR = 0 | (( MQM / MQMR ) > 0.25 & ( MQM / MQMR ) < 1.75 )" \
-f "PAIREDR = 0 | ( PAIRED > 0.05 & PAIREDR > 0.05 & ( PAIREDR / PAIRED ) < 1.75 & ( PAIREDR / PAIRED ) > 0.25 ) " \
-f "SAF > 0 & SAR > 0 & RPL > 1 & RPR > 1" \
-f "AC > 0" | \
vcfallelicprimitives --keep-info --keep-geno | \
vt decompose_blocksub - -o + | \
vt normalize + -m -r $REF -o + | \
bcftools norm --rm-dup all -Ov | \
vcfclassify - | \
vcffilter \
-v -f "INS | DEL" \
-f "AC > 0" | \
bgzip -c > $OUTDIR/$IND\_whole_genome_filter1.vcf.gz;
tabix $OUTDIR/$IND\_whole_genome_filter1.vcf.gz;

## Re-genotype half-calls
python3 $correct_halfcalls -i $OUTDIR/$IND\_whole_genome_filter1.vcf.gz \
-o $OUTDIR/$IND\_whole_genome_filter1_regeno.vcf;
bgzip $OUTDIR/$IND\_whole_genome_filter1_regeno.vcf;
tabix $OUTDIR/$IND\_whole_genome_filter1_regeno.vcf.gz;

## Mean individual depth
vcftools --gzvcf $OUTDIR/$IND\_whole_genome_filter1_regeno.vcf.gz \
--depth --out $OUTDIR/$IND\_whole_genome;
cat $OUTDIR/$IND\_whole_genome.idepth | tail -n 1 | cut -f 3 > $OUTDIR/$IND\_depth;

## Filter 2
vcftools --gzvcf $OUTDIR/$IND\_whole_genome_filter1_regeno.vcf.gz \
--min-meanDP $(echo "scale=2; $(cat $OUTDIR/$IND\_depth) / 3" | bc) \
--max-meanDP $(echo "scale=2; $(cat $OUTDIR/$IND\_depth) * 2" | bc) \
--recode --recode-INFO-all --stdout | \
bgzip -c > $OUTDIR/$IND\_whole_genome_filter2.vcf.gz;
tabix $OUTDIR/$IND\_whole_genome_filter2.vcf.gz;

### Create consensus sequence
cat $REF | bcftools consensus -H R -s $IND $OUTDIR/$IND\_whole_genome_filter2.vcf.gz \
> $OUTDIR0/$PROJECT1\_consensus_reference.fasta;

### Index consensus sequence
bwa index $OUTDIR0/$PROJECT1\_consensus_reference.fasta;
samtools faidx $OUTDIR0/$PROJECT1\_consensus_reference.fasta;
