#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J gone
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

MIN_SCAFFOLD_SIZE=1000000;
PHASE=2 ### Phase = 0 (pseudohaploids), 1 (known phase), 2 (unknown phase)
cMMbA=3.06 # Autosomal recombination rate ### CentiMorgans per Megabase (if distance is not available in map file). 
cMMbZ=1.16 # Z recombination rate
DIST=2 ### none (0), Haldane correction (1) or Kosambi correction (2)
ZERO=0 ### 0: Remove SNPs with zeroes (1: allow for them)
maxNSNP=100000;
MAC=0;
REPS=1000; ### Number of replicates to run GONE (recommended 40)

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
OUTDIR0=$MAINDIR/working/Simon/$PROJECT;
WORKDIR=$MAINDIR/data/$PROJECT;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A9_define_regions;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
METADATA=$WORKDIR/metadata;

### Load modules
module load bioinfo-tools BEDTools/2.29.2 bcftools/1.17 vcflib/1.0.1 vcftools/0.1.16 plink/1.90b4.9;
GONE=/crex/proj/snic2020-2-25/bin//GONE-accessed-2023-02-07/Linux;

mkdir $OUTDIR0/gone;
OUTDIR0=$OUTDIR0/gone;
mkdir $OUTDIR0/metadata;

# Perform msmc analysis
for DATA in $(echo autosomal homogametic); do

  ### Create folders
  mkdir $OUTDIR0/$DATA;
  OUTDIR=$OUTDIR0/$DATA;
  cd $OUTDIR;

  # Copy files
  cp -r $GONE/PROGRAMMES $OUTDIR;
  cp $GONE/script_GONE.sh $OUTDIR;

  # Set up parameter file
  if [ $DATA == "autosomal" ]; then
    cMMb=$cMMbA;
    GENOMES=$(cat $METADATA/sample_info.txt | tail -n+2 | \
    awk -F'\t' 'BEGIN {SUM=0} {SUM+=2} END {print SUM}');
  elif [ $DATA == "homogametic" ]; then
    cMMb=$cMMbZ;
    GENOMES=$(cat $METADATA/sample_info.txt | tail -n+2 | \
    awk -F'\t' 'BEGIN {SUM=0} {if($3=="Female") {SUM+=1} else {SUM+=2}} END {print SUM-1}'); # Remove uneven numbered female with lowest coverage
  fi;

  MAF=$(echo $MAC $GENOMES | awk '{if($1 > 0) {print ($1+0.5)/$2} else {print 0}}');

  cat $GONE/INPUT_PARAMETERS_FILE | \
  awk -v PHASE=$PHASE -v cMMb=$cMMb -v DIST=$DIST -v ZERO=$ZERO -v maxNSNP=$maxNSNP -v MAF=$MAF -v REPS=$REPS \
  '{if($0~/^PHASE/) print "PHASE="PHASE" ### Phase = 0 (pseudohaploids), 1 (known phase), 2 (unknown phase)"; \
  else if($0~/^cMMb/) print "cMMb="cMMb" ### CentiMorgans per Megabase (if distance is not available in map file).";
  else if($0~/^DIST/) print "DIST="DIST" ### none (0), Haldane correction (1) or Kosambi correction (2)"; \
  else if($0~/^MAF/) print "MAF="MAF" ### Minor allele frequency (0-1) (recommended 0)"; \
  else if($0~/^ZERO/) print "ZERO="ZERO" ### 0: Remove SNPs with zeroes (1: allow for them)"; \
  else if($0~/^maxNSNP/) print "maxNSNP="maxNSNP" ### Maximum approx number of SNPs per chromosomes to be analysed (maximum number is 50000)"; \
  else if($0~/^REPS/) print "REPS="REPS" ### Number of replicates to run GONE (recommended 40)"; \
  else print}' \
  > $OUTDIR/INPUT_PARAMETERS_FILE;

  ### Find callable regions (maximum 200 scaffolds)
  for SCAFFOLD in $(cat $BEDS/$PROJECT\_$DATA.bed | cut -f 1 | sort -u); do
  if [ $(cat $BEDS/$PROJECT\_$DATA.bed | grep $SCAFFOLD | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}') -ge $MIN_SCAFFOLD_SIZE ]; then
      cat $BEDS/$PROJECT\_$DATA.bed | grep $SCAFFOLD;
  fi;
  done | awk -F'\t' 'BEGIN{chr="null"; chrnr="0"} {if(chr!=$1) {chr=$1; chrnr++}; if(chrnr < 200) print}' \
  > $OUTDIR0/metadata/$PROJECT\_$DATA.bed;

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
    # Remove uneven numbered female with lowest coverage
    bcftools view $VCFS/$PROJECT\_$DATA.vcf.gz -R $OUTDIR0/metadata/$PROJECT\_$DATA.bed -s^95694 | \
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

  plink --vcf $OUTDIR/$PROJECT\_$DATA.vcf.gz \
  --vcf-half-call missing --allow-extra-chr --recode --out $OUTDIR/$PROJECT\_gone_data;
  mv $OUTDIR/$PROJECT\_gone_data.map $OUTDIR/$PROJECT\_gone_data.map2;
  cat $OUTDIR/$PROJECT\_gone_data.map2 | \
  awk -F'\t' 'BEGIN{chr="null"; chrnr="0"} {if(chr!=$1) {chr=$1; chrnr++} print chrnr" SNP"NR" 0 "$4}' \
  > $OUTDIR/$PROJECT\_gone_data.map;

  # Remove temporary files
  rm $OUTDIR/$PROJECT\_gone_data.map2;

  # run GONE
  bash script_GONE.sh $PROJECT\_gone_data;

done;
