#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J PAR_dxy
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

SCAFFOLD=CADDXX010000137.1;
WINDOW=2000;
STEP=500;
FIRST_RUN=Yes; # If "Yes", then all files will be created from scratch. Otherwise previous files will be reused, and only window-settigns will be updated. Accepts Yes or No

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
module load bioinfo-tools samtools/1.14 python/3.9.5 bcftools/1.17 vcflib/1.0.1 BEDTools/2.29.2 R_packages/4.0.0;
msa2vcf=/crex/proj/snic2020-2-25/bin/msa2vcf-accessed-2023-11-14/msa2vcf.py;
diploid2haploid=/crex/proj/snic2020-2-25/nobackup/simon/bin_tools/diploid2haploid.py;
genomics_general=/crex/proj/snic2020-2-25/bin/genomics_general-accessed-2021-12-09;

# Define functions
genome_scans_windows=$FUNCTIONS/genome_scans_windows_PAR.r;

### Create folders
mkdir $OUTDIR/PAR_dxy;
OUTDIR=$OUTDIR/PAR_dxy;
mkdir $OUTDIR/vcfs \
$OUTDIR/sequences \
$OUTDIR/absolute_divergence;
OUTDIR1=$OUTDIR/absolute_divergence;
mkdir $OUTDIR1/window_$WINDOW\_step_$STEP;

if [ $(echo $FIRST_RUN) == "Yes" ]; then

  ### Set up vcf-files

  ## Full scaffold
  bcftools view $WORKDIR1/A8_PhaseWY/final_output/vcfs/$PROJECT1\_phased_all_variants.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD.vcf.gz;

  bcftools view $WORKDIR2/A8_PhaseWY/final_output/vcfs/$PROJECT2\_phased_all_variants.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD.vcf.gz;

  ## Sex linked scaffold Z
  bcftools view $WORKDIR1/A8_PhaseWY/final_output/vcfs/$PROJECT1\_homogametic_heterogametes.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz;

  bcftools view $WORKDIR2/A8_PhaseWY/final_output/vcfs/$PROJECT2\_homogametic_heterogametes.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz;

  ## Sex linked scaffold W
  bcftools view $WORKDIR1/A8_PhaseWY/final_output/vcfs/$PROJECT1\_heterogametic.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_heterogametic.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_heterogametic.vcf.gz;

  bcftools view $WORKDIR2/A8_PhaseWY/final_output/vcfs/$PROJECT2\_heterogametic.vcf.gz -r $SCAFFOLD | \
  vcfclassify - |  vcfcreatemulti | vcffixup - | \
  vcffilter -f "AC > 0" | vcffilter -s -f "!( INS | DEL )" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_heterogametic.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_heterogametic.vcf.gz;

  ## Convert Project1 to Project2 reference

  # Extract region from species references
  cat $REF1.fai | grep $SCAFFOLD | awk -F'\t' '{print $1":1-"$2}' \
  > $OUTDIR/sequences/$SCAFFOLD\_region.txt;
  samtools faidx $REF1 -r $OUTDIR/sequences/$SCAFFOLD\_region.txt \
  > $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta;
  samtools faidx $REF2 -r $OUTDIR/sequences/$SCAFFOLD\_region.txt \
  > $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta;
  echo -e ">${PROJECT2}" \
  > $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
  samtools faidx $REF2 -r $OUTDIR/sequences/$SCAFFOLD\_region.txt | \
  grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
  >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;

  # Add full scaffold for males
  for IND in $(cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1}'); do

    # Create haplotype sequences for full scaffold
    echo -e ">${IND}_1" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta -H 1pIu $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    echo -e ">${IND}_2" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta -H 2pIu $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;

  done;
  for IND in $(cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1}'); do

    # Create haplotype sequences for full scaffold
    echo -e ">${IND}_1" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta -H 1pIu $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    echo -e ">${IND}_2" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta -H 2pIu $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;

  done;

  # Add sex_linked region of scaffold for females
  for IND in $(cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1}'); do

    # Create haplotype sequences for sex_linked region of scaffold
    echo -e ">${IND}_Z" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta -H 1 $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    echo -e ">${IND}_W" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT1\_$SCAFFOLD\_region.fasta -H 1 $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_heterogametic.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;

  done;
  for IND in $(cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1}'); do

    # Create haplotype sequences for sex_linked region of scaffold
    echo -e ">${IND}_Z" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta -H 1 $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_homogametic_heterogametes.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    echo -e ">${IND}_W" \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;
    samtools faidx $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta | \
    bcftools consensus -s $IND -f $OUTDIR/sequences/$PROJECT2\_$SCAFFOLD\_region.fasta -H 1 $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_heterogametic.vcf.gz | \
    grep -v "^>" |  tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta;

  done;

  # Convert alignment to vcf
  python $msa2vcf $OUTDIR/sequences/$SCAFFOLD\_region_samples.fasta | \
  bcftools view -s ^$PROJECT2 | \
  awk -F'\t' -v SCAFFOLD=$SCAFFOLD \
  '{OFS="\t"; if($0 ~ /^#/) print; \
  else {$1=SCAFFOLD; print}}' | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD.vcf.gz;

  ## Convert genotypes to haploid
  python $diploid2haploid \
  -i $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD.vcf.gz \
  -o $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid.vcf;
  bgzip $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid.vcf;
  tabix $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid.vcf.gz;

  # Fix vcf
  bcftools view $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid.vcf.gz | \
  vcffixup - | \
  vcffilter -f "AC > 0 & AF < 1" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz;
  tabix $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz;

  # Prepare input data
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") {print $1"_1\t1\n"$1"_2\t1"} else {print $1"_Z\t1\n"$1"_W\t1"}}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_all.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") {print $1"_1\t1\n"$1"_2\t1"} else {print $1"_Z\t1\n"$1"_W\t1"}}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_all.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT1=$PROJECT1 '{print $1"_1\t"PROJECT1"\n"$1"_2\t"PROJECT1}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_population_all.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT2=$PROJECT2 '{print $1"_1\t"PROJECT2"\n"$1"_2\t"PROJECT2}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_population_all.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1"_1\t1\n"$1"_2\t1"}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_males.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Male") print $1"_1\t1\n"$1"_2\t1"}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_males.txt;
  cat $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_males.txt | cut -f1 \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_males.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_Z\t1\n"$1"_W\t1"}' \
  > $OUTDIR1/$PROJECT1\_ploidy_females.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_Z\t1\n"$1"_W\t1"}' \
  > $OUTDIR1/$PROJECT2\_ploidy_females.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT1=$PROJECT1 '{if($3=="Male") print $1"_1\t"PROJECT1"\n"$1"_2\t"PROJECT1}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_population_males.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT2=$PROJECT2 '{if($3=="Male") print $1"_1\t"PROJECT2"\n"$1"_2\t"PROJECT2}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_population_males.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_Z\tZ\n"$1"_W\tW"}' \
  > $OUTDIR1/$PROJECT1\_population_females.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_Z\tZ\n"$1"_W\tW"}' \
  > $OUTDIR1/$PROJECT2\_population_females.txt;
  cat $OUTDIR1/$PROJECT1\_ploidy_females.txt | cut -f 1 \
  > $OUTDIR1/$PROJECT1\_females.txt;
  cat $OUTDIR1/$PROJECT2\_ploidy_females.txt | cut -f 1 \
  > $OUTDIR1/$PROJECT2\_females.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_W\t1"}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_females_W.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' '{if($3=="Female") print $1"_W\t1"}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_females_W.txt;
  cat $METADATA1/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT1=$PROJECT1 '{if($3=="Female") print $1"_W\t"PROJECT1}' \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_population_females_W.txt;
  cat $METADATA2/sample_info.txt | tail -n+2 | awk -F'\t' -v PROJECT2=$PROJECT2 '{if($3=="Female") print $1"_W\t"PROJECT2}' \
  >> $OUTDIR1/$PROJECT1\_$PROJECT2\_population_females_W.txt;
  cat $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_females_W.txt | cut -f 1 \
  > $OUTDIR1/$PROJECT1\_$PROJECT2\_females_W.txt;

  python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz \
  --ploidyFile $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_all.txt | \
  bgzip -c > $OUTDIR1/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_snps_all.geno.gz;
  bcftools view $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz \
  -S $OUTDIR1/$PROJECT1\_$PROJECT2\_males.txt | \
  vcffixup - | \
  vcffilter -f "AC > 0 & AF < 1" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_males.vcf.gz;
  python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_males.vcf.gz \
  --ploidyFile $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_males.txt | \
  bgzip -c > $OUTDIR1/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_snps_males.geno.gz;
  bcftools view $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz \
  -S $OUTDIR1/$PROJECT1\_females.txt | \
  vcffixup - | \
  vcffilter -f "AC > 0 & AF < 1" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_haploid_snps_females.vcf.gz;
  python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR/vcfs/$PROJECT1\_$SCAFFOLD\_haploid_snps_females.vcf.gz \
  --ploidyFile $OUTDIR1/$PROJECT1\_ploidy_females.txt | \
  bgzip -c > $OUTDIR1/$PROJECT1\_$SCAFFOLD\_snps_females.geno.gz;
  bcftools view $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz \
  -S $OUTDIR1/$PROJECT2\_females.txt | \
  vcffixup - | \
  vcffilter -f "AC > 0 & AF < 1" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_haploid_snps_females.vcf.gz;
  python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR/vcfs/$PROJECT2\_$SCAFFOLD\_haploid_snps_females.vcf.gz \
  --ploidyFile $OUTDIR1/$PROJECT2\_ploidy_females.txt | \
  bgzip -c > $OUTDIR1/$PROJECT2\_$SCAFFOLD\_snps_females.geno.gz;
  bcftools view $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_all.vcf.gz \
  -S $OUTDIR1/$PROJECT1\_$PROJECT2\_females_W.txt | \
  vcffixup - | \
  vcffilter -f "AC > 0 & AF < 1" | \
  bgzip -c > $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_females_W.vcf.gz;
  python3 $genomics_general/VCF_processing/parseVCF.py -i $OUTDIR/vcfs/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_haploid_snps_females_W.vcf.gz \
  --ploidyFile $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_females_W.txt | \
  bgzip -c > $OUTDIR1/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_snps_females_W.geno.gz;

fi;

### Calculate windows with consideration to missing data
bedtools intersect -a $WORKDIR1/A8_PhaseWY/final_output/beds/$PROJECT1\_target_region.bed \
-b $WORKDIR2/A8_PhaseWY/final_output/beds/$PROJECT2\_target_region.bed | grep $SCAFFOLD \
> $OUTDIR1/$PROJECT1\_$PROJECT2\_callable_all.bed;
Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=0 \
PROJECT=$PROJECT1\_$PROJECT2 REF=$REF2 OUTDIR=$OUTDIR1 DATA=all;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_all.txt | tail -n+2 | cut -f1,2,3 \
> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_all_infile.txt;
bedtools intersect -a $WORKDIR1/A8_PhaseWY/final_output/beds/$PROJECT1\_heterogametic.bed \
-b $WORKDIR2/A8_PhaseWY/final_output/beds/$PROJECT2\_heterogametic.bed | grep $SCAFFOLD \
> $OUTDIR1/$PROJECT1\_$PROJECT2\_callable_heterogametic.bed;
Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=0 \
PROJECT=$PROJECT1\_$PROJECT2 REF=$REF2 OUTDIR=$OUTDIR1 DATA=heterogametic;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic.txt | tail -n+2 | cut -f1,2,3 \
> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;
cat $WORKDIR1/A8_PhaseWY/final_output/beds/$PROJECT1\_heterogametic.bed  | grep $SCAFFOLD \
> $OUTDIR1/$PROJECT1\_callable_heterogametic.bed;
Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=0 \
PROJECT=$PROJECT1 REF=$REF2 OUTDIR=$OUTDIR1 DATA=heterogametic;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_window_$WINDOW\_step_$STEP\_heterogametic.txt | tail -n+2 | cut -f1,2,3 \
> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;
cat $WORKDIR2/A8_PhaseWY/final_output/beds/$PROJECT2\_heterogametic.bed  | grep $SCAFFOLD \
> $OUTDIR1/$PROJECT2\_callable_heterogametic.bed;
Rscript $genome_scans_windows --args WINDOW=$WINDOW STEP=$STEP MIN_SCAFFOLD_SIZE=0 \
PROJECT=$PROJECT2 REF=$REF2 OUTDIR=$OUTDIR1 DATA=heterogametic;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic.txt | tail -n+2 | cut -f1,2,3 \
> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;

## Caluclate pairwise dxy males Z & PAR
python3 $genomics_general/popgenWindows.py -g $OUTDIR1/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_snps_males.geno.gz \
-f phased --ploidyFile $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_males.txt \
-o $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_pairwise_males_window_$WINDOW\_step_$STEP.txt \
--popsFile $OUTDIR1/$PROJECT1\_$PROJECT2\_population_males.txt \
-p $PROJECT1 -p $PROJECT2 --windType predefined --writeFailedWindows \
--windCoords $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_all_infile.txt;

## Caluclate pairwise dxy Z W project 1
python3 $genomics_general/popgenWindows.py -g $OUTDIR1/$PROJECT1\_$SCAFFOLD\_snps_females.geno.gz \
-f phased --ploidyFile $OUTDIR1/$PROJECT1\_ploidy_females.txt \
-o $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$SCAFFOLD\_pairwise_females_window_$WINDOW\_step_$STEP.txt \
--popsFile $OUTDIR1/$PROJECT1\_population_females.txt \
-p Z -p W --windType predefined --writeFailedWindows \
--windCoords $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;

## Caluclate pairwise dxy Z W project 2
python3 $genomics_general/popgenWindows.py -g $OUTDIR1/$PROJECT2\_$SCAFFOLD\_snps_females.geno.gz \
-f phased --ploidyFile $OUTDIR1/$PROJECT2\_ploidy_females.txt \
-o $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT2\_$SCAFFOLD\_pairwise_females_window_$WINDOW\_step_$STEP.txt \
--popsFile $OUTDIR1/$PROJECT2\_population_females.txt \
-p Z -p W --windType predefined --writeFailedWindows \
--windCoords $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;

## Caluclate pairwise dxy female W
python3 $genomics_general/popgenWindows.py -g $OUTDIR1/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_snps_females_W.geno.gz \
-f phased --ploidyFile $OUTDIR1/$PROJECT1\_$PROJECT2\_ploidy_females_W.txt \
-o $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_pairwise_females_W_window_$WINDOW\_step_$STEP.txt \
--popsFile $OUTDIR1/$PROJECT1\_$PROJECT2\_population_females_W.txt \
-p $PROJECT1 -p $PROJECT2 --windType predefined --writeFailedWindows \
--windCoords $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_window_$WINDOW\_step_$STEP\_heterogametic_infile.txt;

## Combine data
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_pairwise_males_window_$WINDOW\_step_$STEP.txt | \
tr ',' '\t' | awk -F="\t" '{if(NR==1) {print $0"\tdata_type"} else {print $0"\tmales"}}' \
> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_dxy_all_data_window_$WINDOW\_step_$STEP.txt;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$SCAFFOLD\_pairwise_females_window_$WINDOW\_step_$STEP.txt | \
tr ',' '\t' | awk -F="\t" -v PROJECT1=$PROJECT1 '{print $0"\tfemales_"PROJECT1}' | tail -n+2 \
>> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_dxy_all_data_window_$WINDOW\_step_$STEP.txt;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT2\_$SCAFFOLD\_pairwise_females_window_$WINDOW\_step_$STEP.txt | \
tr ',' '\t' | awk -F="\t" -v PROJECT2=$PROJECT2 '{print $0"\tfemales_"PROJECT2}' | tail -n+2 \
>> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_dxy_all_data_window_$WINDOW\_step_$STEP.txt;
cat $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_pairwise_females_W_window_$WINDOW\_step_$STEP.txt | \
tr ',' '\t' | awk -F="\t" '{print $0"\tfemales_W"}' | tail -n+2 \
>> $OUTDIR1/window_$WINDOW\_step_$STEP/$PROJECT1\_$PROJECT2\_$SCAFFOLD\_dxy_all_data_window_$WINDOW\_step_$STEP.txt;
