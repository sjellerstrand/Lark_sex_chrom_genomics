#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J linkage_decay
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

MIN_SCAFFOLD_SIZE=1000000;
SNPS_DENSITY=1; ### SNPs/kbp
BINS=1000; # Number of bp to average bins over
MAC=1;

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
PROJECT2=Skylark_2021;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT;
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A9_define_regions;
METADATA=$WORKDIR/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;

### Load modules
module load bioinfo-tools bcftools/1.17 vcflib/1.0.1 vcftools/0.1.16 BEDTools/2.29.2 python/3.9.5 R_packages/4.0.0;
PopLDdecay=/proj/snic2020-2-25/bin/PopLDdecay-accessed-2023-07-16/bin;

# Define functions
LD_decay_windows=$FUNCTIONS/LD_decay_windows_PAR.r

### Create folders
mkdir $OUTDIR0/E1_linkage_decay;
OUTDIR0=$OUTDIR0/E1_linkage_decay;
mkdir $OUTDIR0/metadata;

# Calculate Linkage decay
for DATA in $(echo autosomal homogametic); do

  ### Create folders
  mkdir $OUTDIR0/$DATA;
  OUTDIR=$OUTDIR0/$DATA;
  cd $OUTDIR;

  # Set up parameter file
  if [ $DATA == "autosomal" ]; then
    GENOMES=$(cat $METADATA/sample_info.txt | tail -n+2 | \
    awk -F'\t' 'BEGIN {SUM=0} {SUM+=2} END {print SUM}');
  elif [ $DATA == "homogametic" ]; then
    GENOMES=$(cat $METADATA/sample_info.txt | tail -n+2 | \
    awk -F'\t' 'BEGIN {SUM=0} {if($3=="Female") {SUM+=1} else {SUM+=2}} END {print SUM}');
  fi;

  ### Find callable regions
  for SCAFFOLD in $(cat $BEDS/$PROJECT\_$DATA.bed | cut -f 1 | sort -u); do
  if [ $(cat $REF.fai | grep $SCAFFOLD | cut -f2) -ge $MIN_SCAFFOLD_SIZE ]; then
      cat $BEDS/$PROJECT\_$DATA.bed | grep $SCAFFOLD;
  fi;
  done > $OUTDIR0/metadata/$PROJECT\_$DATA.bed;

  # Set up input files
  if [ $DATA == "autosomal" ]; then
    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz -R $OUTDIR0/metadata/$PROJECT\_$DATA.bed | \
    vcfclassify - | \
    vcffixup - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > $MAC & AC < ( $GENOMES - 1 )" | \
    vcftools --gzvcf - --max-alleles 2 --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA.vcf.gz;

  elif [ $DATA == "homogametic" ]; then

    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz -R $OUTDIR0/metadata/$PROJECT\_$DATA.bed | \
    vcfclassify - | \
    vcffixup - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > $MAC & AC < ( $GENOMES - 1 )" | \
    vcftools --gzvcf - --max-alleles 2 --recode --recode-INFO-all --stdout | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA\_temp.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA\_temp.vcf.gz;

    bcftools view $OUTDIR/$PROJECT\_$DATA\_temp.vcf.gz \
    -s $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}' | tr '\n' ',' | sed 's/.$//') | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA\_males_temp.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA\_males_temp.vcf.gz;

    bcftools view $OUTDIR/$PROJECT\_$DATA\_temp.vcf.gz \
    -s $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}' | grep -v "95694" | tr '\n' ',' | sed 's/.$//') | \
    awk 'BEGIN {OFS="\t"} {if($0 ~ /^#CHROM/) {
          for(i=10; i<=NF; i+=2) {$i=$i$(i+1); $(i+1)=""} print $0}
        else if($0 ~ /^#/) {print}
        else {for(i=10; i<=NF; i+=2) {$i=$i"|"$(i+1); $(i+1)=""} print $0}}' | \
    tr -s '\t' | sed 's/\t$//' | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA\_females_temp.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA\_females_temp.vcf.gz;

    bcftools merge $OUTDIR/$PROJECT\_$DATA\_males_temp.vcf.gz \
    $OUTDIR/$PROJECT\_$DATA\_females_temp.vcf.gz | \
    vcffixup - | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA.vcf.gz;

    rm $OUTDIR/*temp.vcf.gz*;

  fi;

  CALLABLE=$(cat $OUTDIR0/metadata/$PROJECT\_$DATA.bed | awk -F'\t' 'BEGIN {SUM=0} {SUM+=$3-$2} END {print SUM}');
  NUM_SNPS=$(echo "scale=0; ( $CALLABLE * $SNPS_DENSITY ) / 1000" | bc);
  bcftools view $OUTDIR/$PROJECT\_$DATA.vcf.gz -H | \
  awk -F'\t' '{print $1"\t"$2-1"\t"$2}' | \
  shuf -n $NUM_SNPS \
  > $OUTDIR0/metadata/$PROJECT\_$DATA\_target_pos_subsample.bed;

  bcftools view $OUTDIR/$PROJECT\_$DATA.vcf.gz \
  -R $OUTDIR0/metadata/$PROJECT\_$DATA\_target_pos_subsample.bed | \
  bgzip -c > $OUTDIR/$PROJECT\_$DATA\_sub.vcf.gz;
  tabix $OUTDIR/$PROJECT\_$DATA\_sub.vcf.gz;

  $PopLDdecay/PopLDdecay -InVCF $OUTDIR/$PROJECT\_$DATA\_sub.vcf.gz \
  -MaxDist 10000 -MAF 0 -Het 1 -Miss 1 \
  -OutStat $OUTDIR/$PROJECT\_$DATA.LDdecay -OutType 1;

  # Average over bins
  zcat $OUTDIR/$SEX/$PROJECT\_$DATA.LDdecay.stat.gz | \
  awk -F'\t' -v BINS=$BINS 'BEGIN {SUM_r=0; Sum_p=0} {if(NR==1) {print "Min\tMax\tMid\tMean_r2"} \
  else {SUM_r=+$4; SUM_p=+$6; if((NR-1)%BINS==0)  {print $1-BINS+1"\t"$1"\t"$1-((BINS-1)/2)"\t"SUM_r/SUM_p}}}' \
  > $OUTDIR/$PROJECT\_$DATA.LDdecay.stat.bins.tsv;

done;
