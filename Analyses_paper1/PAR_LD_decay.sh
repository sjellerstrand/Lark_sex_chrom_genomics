#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J PAR_LD_decay
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

SCAFFOLD=CADDXX010000137.1;
MAX_SNPS=2000;

# Parameters to vary:
FIRST_RUN=Yes; # If "Yes", then all files will be created from scratch. Otherwise previous files will be reused, and only window-settigns will be updated. Accepts Yes or No
WINDOW=5000;
STEP=1250;

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
PROJECT1=Rasolark_2021;
PROJECT2=Skylark_2021;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
REF1=$WORKDIR1/A4_make_new_reference/$PROJECT1\_consensus_reference.fasta;
REF2=$WORKDIR2/A4_make_new_reference/$PROJECT2\_consensus_reference.fasta;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
FUNCTIONS=$MAINDIR/scripts/$PROJECT/analyses;

### Load modules
module load bioinfo-tools bcftools/1.17 vcflib/1.0.1 BEDTools/2.29.2 python/3.9.5 R_packages/4.0.0;
PopLDdecay=/proj/snic2020-2-25/bin/PopLDdecay-accessed-2023-07-16/bin;

# Define functions
LD_decay_windows=$FUNCTIONS/LD_decay_windows_PAR.r;

### Create folders
mkdir $OUTDIR/PAR_LD_decay;
OUTDIR=$OUTDIR/PAR_LD_decay;
mkdir $OUTDIR/vcfs \
$OUTDIR/LD_decay;
OUTDIR1=$OUTDIR/LD_decay;
mkdir $OUTDIR1/window_$WINDOW\_step_$STEP;

if [ $(echo $FIRST_RUN) == "Yes" ]; then

  # Set up sample info
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1}' \
  > $OUTDIR/$PROJECT1\_males.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1}' \
  > $OUTDIR/$PROJECT2\_males.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1}' \
  > $OUTDIR/$PROJECT1\_females.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1}' \
  > $OUTDIR/$PROJECT2\_females.txt;

  ## Set up vcf-files

  for PROJ in $(echo 1 2); do
    for SEX in $(echo males females); do

      WORKDIR3=$(eval "echo \${WORKDIR$PROJ}");
      PROJECT3=$(eval "echo \${PROJECT$PROJ}");

      # Extract full scaffold
      bcftools view $WORKDIR3/A8_PhaseWY/final_output/vcfs/$PROJECT3\_phased_all_variants.vcf.gz \
      -r $SCAFFOLD -S $OUTDIR/$PROJECT3\_$SEX.txt | \
      vcfclassify - | vcfcreatemulti | vcffixup - | \
      vcffilter -f "AC > 0 & AF < 1" | vcffilter -s -f "!( INS | DEL )" | \
      awk -F'\t' '{OFS="\t"; if($0 ~ /^#/) {print} else {if($5 !~ /,/) print}}' | \
      awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
      awk -F'\t' '{OFS="\t"; if($0 ~ /^#/) {print} else {$3=$1"_"$2; print}}' | \
      bgzip -c > $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz;
      tabix $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz;

      # Find longest non-variable distance
      bcftools view -H $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz | \
      cut -f 2 | awk '{if($NR==1) {i=$1} else {print $1-i-1; i=$1}}' | sort -nr | head -n1 \
      >>  $OUTDIR/vcfs/max_non_var_dist.txt;

      # Find number of snps
      bcftools view $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz -H | wc -l \
      >> $OUTDIR/vcfs/minimum_snps.txt;

    done;
  done;

  ### Downsample variants to uniform disribution

  ## Set bins and number of snps per bin
  cat $REF1.fai | grep $SCAFFOLD \
  > $OUTDIR/vcfs/$SCAFFOLD.fai;
  SCAFFOLD_LENGTH=$(cat $REF1.fai | grep $SCAFFOLD | cut -f 2);
  WIND_SIZE=$(cat $OUTDIR/vcfs/max_non_var_dist.txt | sort -n | head -n1);
  NB_SNPS=$(cat $OUTDIR/vcfs/minimum_snps.txt | sort -n | head -n1);
  WIND_SNPS=$(echo $NB_SNPS $SCAFFOLD_LENGTH $WIND_SIZE | awk '{print int(($1/($2/$3))*0.5 + 0.5)}');
  WIND_SNPS=$( echo $WIND_SNPS $SCAFFOLD_LENGTH $WIND_SIZE $MAX_SNPS | \
  awk '{if(int($1*($2/$3)) < $4) {print $1} else {print int($4/($2/$3) + 0.5)}}');
  bedtools makewindows -g $OUTDIR/vcfs/$SCAFFOLD.fai -w $WIND_SIZE | \
  awk -F'\t' '{print $1":"$2+1"-"$3}' \
  > $OUTDIR/vcfs/$SCAFFOLD\_windows.txt;

  ## Downsample scaffold
  for PROJ in $(echo 1 2); do
    for SEX in $(echo males females); do

      WORKDIR3=$(eval "echo \${WORKDIR$PROJ}");
      PROJECT3=$(eval "echo \${PROJECT$PROJ}");

      ## Extract variants per window
      bcftools view $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz -h \
      > $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform.vcf;

      for i in $(seq 1 1 $(cat $OUTDIR/vcfs/$SCAFFOLD\_windows.txt | wc -l)); do
        bcftools view $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps.vcf.gz -H \
        -r $(cat $OUTDIR/vcfs/$SCAFFOLD\_windows.txt | head -n$i | tail -n1) | \
        shuf | head -n$WIND_SNPS \
        >> $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform.vcf;
      done;

      ## Sort and subsample
      bcftools sort $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform.vcf | \
      bgzip -c > $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform_downsampled.vcf.gz;
      tabix $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform_downsampled.vcf.gz;
      rm $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform.vcf;

      # Calculate Linkage decay
      $PopLDdecay/PopLDdecay -InVCF $OUTDIR/vcfs/$PROJECT3\_$SCAFFOLD\_$SEX\_snps_uniform_downsampled.vcf.gz \
      -MaxDist 100 -MAF 0 -Het 1 -Miss 1 \
      -OutStat $OUTDIR1/$PROJECT3\_$SCAFFOLD\_$SEX.LDdecay -OutType 8;

    done
  done;

fi;

## Calcualte windowed averages
for PROJ in $(echo 1 2); do
  for SEX in $(echo males females); do

    WORKDIR3=$(eval "echo \${WORKDIR$PROJ}");
    PROJECT3=$(eval "echo \${PROJECT$PROJ}");

    # Calculate windowed means
    Rscript $LD_decay_windows --args WINDOW=$WINDOW STEP=$STEP OUTDIR=$OUTDIR1 PROJECT=$PROJECT3 DATA=$SEX \
    SCAFFOLD=$SCAFFOLD SCAFFOLD_LENGTH=$SCAFFOLD_LENGTH;

    ## Combine data
    cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT3\_$SCAFFOLD\_$SEX\_window_$WINDOW\_step_$STEP.txt | \
    awk -F="\t" -v PROJECT=$PROJECT3 -v SEX=$SEX '{if(NR==1) {print $0"\tdata_type"} else {print $0"\t"PROJECT"_"SEX}}' \
    >> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_LD_decay_all_data_window_$WINDOW\_step_$STEP.txt;

  done;
done;

