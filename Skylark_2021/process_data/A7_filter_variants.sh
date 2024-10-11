#!/bin/bash -l

#SBATCH -A naiss2023-5-169
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 150:00:00
#SBATCH -J filter_variants
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

SEX_SYSTEM=ZW;
NONDIPLOIDS=No;
EXCESS_HETEROZYGOSITY=Yes; ### Filter sites that are heterozygos in all diploid individuals

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
PROJECT2=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
VCF_IN=$WORKDIR/A6_call_variants/$PROJECT\_raw;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT;

### Load modules
module load bioinfo-tools vcflib/1.0.1 vt/0.5772 bcftools/1.14 \
vcftools/0.1.16 python/3.9.5 R_packages/4.0.0 plink/1.90b4.9;

### Define functions
filter_stats=$FUNCTIONS/quality_control/filter_stats.sh;
pca=$FUNCTIONS/quality_control/pca.sh;
contam=$FUNCTIONS/quality_control/contam.sh;

### Create folders
mkdir $OUTDIR0/A7_filter_variants;
OUTDIR0=$OUTDIR0/A7_filter_variants;

### Apply filters

## Filter 0: Quality check raw vcf
filter=0;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/$PROJECT\_$filter;
cp $VCF_IN.vcf.gz $VCF_OUT.vcf.gz;
cp $VCF_IN.vcf.gz.tbi $VCF_OUT.vcf.gz.tbi;
source $filter_stats;
source $pca;
source $contam;

## Filter 1: Hard filters, individual alleles, and re-genotyping
VCF_IN=$VCF_OUT;
filter=1;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/$PROJECT\_$filter;
vcftools --gzvcf $VCF_IN.vcf.gz \
--exclude-bed $WORKDIR2/B1_prepare_reference/$PROJECT2\_masked_repeats.bed \
--max-missing 0.8 \
--min-meanDP 5 \
--max-meanDP $(echo "scale=2; $(cat $VCF_IN\_stats.txt | grep "var mean" | cut -f4) * 3" | bc) \
--minQ 30 \
--minDP 5 \
--recode --recode-INFO-all --stdout | \
vcffilter \
-f "QA = 30 | QA > 30" \
-f "MQMR = 0 | (( MQM / MQMR ) > 0.25 & ( MQM / MQMR ) < 1.75 )" \
-f "PAIREDR = 0 | ( PAIRED > 0.05 & PAIREDR > 0.05 & ( PAIREDR / PAIRED ) < 1.75 & ( PAIREDR / PAIRED ) > 0.25 ) " \
-f "SAF > 0 & SAR > 0 & RPL > 1 & RPR > 1" \
-f "AC > 0" | \
bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $contam;

## Filter 2: Final depth and missingness + excess heterozygosity
VCF_IN=$VCF_OUT;
filter=2;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/$PROJECT\_$filter;

# Get individual mean depths
cat $METADATA/sample_info.txt | tail -n+2 | cut -f1,3 | sort -k1> $OUTDIR/INDS1.txt;
cat $VCF_IN.idepth | tail -n+2 | cut -f1,3 | sort -k1 > $OUTDIR/INDS2.txt;
join $OUTDIR/INDS1.txt $OUTDIR/INDS2.txt | tr ' ' '\t' > $OUTDIR/INDS.txt;

# Get min mean depth
if [ $SEX_SYSTEM != NA ]; then
     if [ $SEX_SYSTEM == ZW ]; then
         cat $OUTDIR/INDS.txt | awk -F'\t' '{if($2=="Female") print $3/2; else print $3}' \
         > $OUTDIR/CORRECTED_INDS_DEPTH.txt;
     elif [ $SEX_SYSTEM == XY ]; then
         cat $OUTDIR/INDS.txt | awk -F'\t' '{if($2=="Male") print $3/2; else print $3}' \
         > $OUTDIR/CORRECTED_INDS_DEPTH.txt;
     fi;
     MIN_MEAN=$(cat $OUTDIR/CORRECTED_INDS_DEPTH.txt | awk 'BEGIN {sum=0; n=0} {sum+=$1; n+=1} END {print sum/(3*n)}');
else
     MIN_MEAN=$(cat $OUTDIR/INDS.txt | awk -F'\t' 'BEGIN {sum=0; n=0} {sum+=$3; n+=1} END {print sum/(3*n)}');
fi;
if [ $(echo "$MIN_MEAN < 5" | bc) -eq 1 ]; then
    MIN_MEAN=5;
fi;
echo "Applying min_meanDP =" $MIN_MEAN;

# Get max mean depth
MAX_MEAN=$(echo "scale=2; $(cat $VCF_IN\_stats.txt | grep "var mean" | cut -f4) * 1.5" | bc)
echo "Applying max-meanDP =" $MAX_MEAN;

# Get min depth
if [ $SEX_SYSTEM != NA ]; then
     MIN_DP=$(echo "scale=2; $(cat $OUTDIR/CORRECTED_INDS_DEPTH.txt | sort -k1 -n | head -n1) / 3" | bc);
else
     MIN_DP=$(echo "scale=2; $(cat $OUTDIR/INDS.txt | cut -f 3 | sort -k1 -n | head -n1) / 3" | bc);
fi;
if [ $(echo "$MIN_DP < 5" | bc) -eq 1 ]; then
    MIN_DP=5;
fi;
echo "Applying minDP =" $MIN_DP;

# Get missingness
MISSING=$(if [ $(echo "$(cat $VCF_IN\_stats.txt | grep "var miss" | cut -f 6) < 0.01" | bc) -eq 1 ]; then echo 0.95; else echo 0.9; fi;);
echo "Applying missing =" $MISSING;

# Save paramters in file
echo -e "MIN_MEAN=$MIN_MEAN\n\
MAX_MEAN=$MAX_MEAN\n\
MIN_DP=$MIN_DP\n\
MISSING=$MISSING"\
> $OUTDIR/filtering_parameters.txt;

# Remove excess heterozygosity
if [ $EXCESS_HETEROZYGOSITY == Yes ]; then
  if [ $NONDIPLOIDS == Yes ]; then
    DIPLOIDS=$(cat $METADATA/ploidy.txt | awk -F'\t' '{if($2==1) print "--remove-indv",$1}' | tr '\n' ' ');
    vcftools --gzvcf $VCF_IN.vcf.gz --hardy $DIPLOIDS --stdout;
  else
    vcftools --gzvcf $VCF_IN.vcf.gz --hardy --stdout;
  fi | tail -n+2 | cut -f1,2,3 | awk -F'\t|/' 'BEGIN {print "#Header"} {if($3 == 0 && $5 == 0) print $1"\t"$2-1"\t"$2}' \
  > $OUTDIR/excess_het_with_header.bed;
else
  echo "#Header" > $OUTDIR/excess_het_with_header.bed;
fi;

# Apply filters and update AC, AF and AN in the INFO field
vcftools --gzvcf $VCF_IN.vcf.gz \
--min-meanDP $MIN_MEAN \
--max-meanDP $MAX_MEAN \
--minDP $MIN_DP \
--exclude-bed $OUTDIR/excess_het_with_header.bed \
--recode --recode-INFO-all --stdout | \
vcffixup - | \
vcffilter \
-f "AC > 0" | \
vcftools --vcf - \
--max-missing $MISSING \
--recode --recode-INFO-all --stdout | \
bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $contam;

## Filter 3: Decompose variants
VCF_IN=$VCF_OUT;
filter=3;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/$PROJECT\_$filter;
vcfallelicprimitives --keep-info --keep-geno $VCF_IN.vcf.gz | \
vt decompose_blocksub - -o + | \
vt normalize + -m -r $REF -o + | \
bcftools norm --rm-dup all -Ov | \
awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
vcffixup - | \
vcfclassify - | \
vcffilter \
-f "AC > 0" | \
vcftools --vcf - \
--max-missing $MISSING \
--exclude-bed $WORKDIR2/B1_prepare_reference/$PROJECT2\_masked_repeats.bed \
--recode --recode-INFO-all --stdout | \
bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $contam;
