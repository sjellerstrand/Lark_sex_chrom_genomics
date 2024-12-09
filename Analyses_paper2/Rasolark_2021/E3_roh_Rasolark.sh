#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J roh
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

MIN_SCAFFOLD_SIZE=1000000;
REC_RATE_A=3.06; # Autosomal recombination rate
REC_RATE_Z=1.16; # Z recombination rate

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT;
WORKDIR=$MAINDIR/data/$PROJECT;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A9_define_regions;
METADATA=$WORKDIR/metadata;

### Load modules
module load bioinfo-tools BEDTools/2.29.2 bcftools/1.17 vcflib/1.0.1 vcftools/0.1.16;

### Create folders
mkdir $OUTDIR0/E3_roh;
OUTDIR0=$OUTDIR0/E3_roh;
mkdir $OUTDIR0/metadata;

# Analyse runs of homozygosity
for DATA in $(echo autosomal homogametic); do

  ### Create folders
  mkdir $OUTDIR0/$DATA;
  OUTDIR=$OUTDIR0/$DATA;

  ### Create mask
  for SCAFFOLD in $(cat $BEDS/$PROJECT\_$DATA.bed | cut -f 1 | sort -u); do
    if [ $(cat $REF.fai | grep $SCAFFOLD | cut -f2) -ge $MIN_SCAFFOLD_SIZE ]; then
      cat $BEDS/$PROJECT\_$DATA.bed | grep $SCAFFOLD;
    fi;
  done > $OUTDIR0/metadata/$PROJECT\_$DATA.bed;
  MASK=$OUTDIR0/metadata/$PROJECT\_$DATA.bed;

  ### Get total region size
  REG_SIZE=$(cat $MASK |awk -F'\t' 'BEGIN {SUM=0} {SUM+=$3-$2} END {print SUM}');

  if [ $DATA == "autosomal" ]; then

    REC_RATE=$REC_RATE_A;

    ### Prepare vcf
    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz -R $MASK | \
    vcfclassify - | \
    vcffixup - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > 0 & AF < 1" | \
    bgzip -c > $OUTDIR/$PROJECT\_$DATA.vcf.gz;
    tabix $OUTDIR/$PROJECT\_$DATA.vcf.gz;

  elif [ $DATA == "homogametic" ]; then 

    REC_RATE=$REC_RATE_Z;

    # Remove uneven numbered female with lowest coverage
    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz -R $OUTDIR0/metadata/$PROJECT\_$DATA.bed -s^95694 | \
    vcfclassify - | \
    vcffixup - | \
    vcffilter -s -f "!( INS | DEL | MNP )" -f "AC > 0 & AF < 1" | \
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

  # Estiamte runs og homozygosity
  bcftools roh -R $MASK $OUTDIR/$PROJECT\_$DATA.vcf.gz \
  -G0 --rec-rate $REC_RATE --AF-tag AF -I -O r | \
  awk -F'\t' '{if($6 >= 50000) print}' | cut -f2- | \
  awk -F'\t' -v REC_RATE=$REC_RATE \
  '{if(NR==1) {print "Sample\tScaffold\tStart\tEnd\tLength\tN_SNPs\tQuality\tAge_generations"} \
  else {AGE=100/(($5/1000000)*2*REC_RATE); print $0"\t"AGE}}' \
  > $OUTDIR/$PROJECT\_$DATA\_roh;

done;
