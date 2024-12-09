#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J snpEff_individuals
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT2;
GENES=$WORKDIR1/D6_snpEff_substitutions/metadata/genes_in.tsv;
FEATURES=$WORKDIR1/D6_snpEff_substitutions/metadata/features_shared.tsv;
DB_NAME=Tgut;
DB_REF=$MAINDIR/data/Skylark_2021/B3_annotation_lift_over/$DB_NAME.fa;
DB_GTF=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.gtf;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR1}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR1}/D1_filter_genes/Falb_uniqID.gtf");

### Load modules
module load bioinfo-tools SeqKit/0.15.0 bcftools/1.17 vcflib/1.0.1 BEDTools/2.29.2 snpEff/4.3t python/3.9.5 gnuparallel/20180822;
msa2vcf=/crex/proj/snic2020-2-25/bin/msa2vcf-accessed-2023-11-14/msa2vcf.py;
diploid2haploid=/crex/proj/snic2020-2-25/nobackup/simon/bin_tools/diploid2haploid.py;
snpeff_config=/sw/bioinfo/snpEff/4.3t/rackham/snpEff.config.orig;

### Activate conda environment
conda activate liftover;

### Create folders
mkdir $OUTDIR/D7_snpEff_individuals;
OUTDIR=$OUTDIR/D7_snpEff_individuals;
mkdir $OUTDIR/metadata \
$OUTDIR/shared_genes \
$OUTDIR/shared_genes/autosomal_autosomal \
$OUTDIR/shared_genes/sex_dropout_sex_dropout \
$OUTDIR/shared_genes/sex_phase_sex_phase \
$OUTDIR/shared_genes/sex_phase_sex_dropout \
$OUTDIR/shared_genes/sex_dropout_sex_phase \
$OUTDIR/shared_genes/sex_phase_autosomal \
$OUTDIR/merge_files \
$OUTDIR/merge_files/autosomal_autosomal \
$OUTDIR/merge_files/sex_dropout_sex_dropout \
$OUTDIR/merge_files/sex_phase_sex_phase \
$OUTDIR/merge_files/sex_phase_sex_dropout \
$OUTDIR/merge_files/sex_dropout_sex_phase \
$OUTDIR/merge_files/sex_phase_autosomal \
$OUTDIR/final_input \
$OUTDIR/data \
$OUTDIR/data/$DB_NAME \
$OUTDIR/data/genomes \
$OUTDIR/$PROJECT1\_autosomal \
$OUTDIR/$PROJECT1\_sexchrom \
$OUTDIR/$PROJECT2\_autosomal \
$OUTDIR/$PROJECT2\_sexchrom \
$OUTDIR/genotype_info \
$OUTDIR/all_data;


# Set up metadata files
OUTGROUPS=$(echo $OUTGROUPS | tr ':' '\n' | cut -d' ' -f1 | tr '\n' '|');
OUTGROUPS=${OUTGROUPS::-1};

# Dropout genes
cat $GENES | awk -F'\t' '{if($3 == "sex_dropout" || $3 == "partial_dropout") print}' | cut -f1 \
> $OUTDIR/metadata/$PROJECT1\_dropout.txt;
cat $GENES | awk -F'\t' '{if($4 == "sex_dropout" || $4 == "partial_dropout") print}' | cut -f1 \
> $OUTDIR/metadata/$PROJECT2\_dropout.txt;

## Per Species files
echo $DB_NAME"|"$OUTGROUPS | tr "|" "\n" \
> $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt;
cp $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt;
cp $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt;
cp $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_A_sampleID.txt;
cp $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_Z_sampleID.txt;
cp $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_W_sampleID.txt;

for IND in $(cat $METADATA1/sample_info.txt | tail -n+2 | cut -f1); do
  echo $IND\_A1;
  echo $IND\_A2;
done >> $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt;
for IND in $(cat $METADATA2/sample_info.txt | tail -n+2 | cut -f1); do
  echo $IND\_A1;
  echo $IND\_A2;
done >> $OUTDIR/metadata/$PROJECT2\_individuals_A_sampleID.txt;
for IND in $(cat $METADATA1/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_Z;
done >> $OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt;
for IND in $(cat $METADATA1/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do
  echo $IND\_Z1;
  echo $IND\_Z2;
done >> $OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt;
for IND in $(cat $METADATA2/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_Z;
done >> $OUTDIR/metadata/$PROJECT2\_individuals_Z_sampleID.txt;
for IND in $(cat $METADATA2/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do
  echo $IND\_Z1;
  echo $IND\_Z2;
done >> $OUTDIR/metadata/$PROJECT2\_individuals_Z_sampleID.txt;
for IND in $(cat $METADATA1/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_W;
done >> $OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt;
for IND in $(cat $METADATA2/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do
  echo $IND\_W;
done >> $OUTDIR/metadata/$PROJECT2\_individuals_W_sampleID.txt;

# Shared species files
cat $OUTDIR/metadata/$PROJECT1\_individuals_A_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_A_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_A_A_sampleID.txt;
cat $OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_Z_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_Z_Z_sampleID.txt;
cat $OUTDIR/metadata/individuals_Z_Z_sampleID.txt \
$OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_W_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_ZW_ZW_sampleID.txt;
cat $OUTDIR/metadata/individuals_Z_Z_sampleID.txt \
$OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_ZW_Z_sampleID.txt;
cat $OUTDIR/metadata/individuals_Z_Z_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_W_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_Z_ZW_sampleID.txt;
cat $OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt \
$OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt \
$OUTDIR/metadata/$PROJECT2\_individuals_A_sampleID.txt | sort -u \
> $OUTDIR/metadata/individuals_ZW_A_sampleID.txt;

export WORKDIR1 WORKDIR2 OUTDIR PROJECT1 PROJECT2 DB_NAME METADATA1 METADATA2 msa2vcf diploid2haploid \
GENES FEATURES OUTGROUPS;

### Convert sequences to database species reference
shared_to_DB_ref() {

  GENE=$1;
  TRANS=$2;
  REGION0_1=$3;
  REGION0_2=$4;

  if [ $REGION0_1 == "sex_phase" ] || [ $REGION0_1 == "partial_dropout" ]; then
    REGION1="sex_phase";
  else
    REGION1=$REGION0_1;
  fi;
  if [ $REGION0_2 == "sex_phase" ] || [ $REGION0_2 == "partial_dropout" ]; then
    REGION2="sex_phase";
  else
    REGION2=$REGION0_2;
  fi;

  # Choose regions
  if [ $REGION1 == "autosomal" ] && [ $REGION2 == "autosomal" ]; then
    SUFFIX="A_A";
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="Z_Z";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="ZW_ZW";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="ZW_Z";
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="Z_ZW";
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "autosomal" ]; then
    SUFFIX="ZW_A";
  fi;

  echo $GENE;
  mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;

  # Find gene CDS features
  cat $FEATURES | grep $GENE\_$TRANS \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv;

  # Sort CDS sequences by exon order
  for EXON in $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv | cut -f17); do
    echo $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv | grep $EXON) $(echo $EXON | rev | cut -d'_' -f2 | rev);
  done | tr ' ' '\t' | sort -nk18,18 | cut -f1-17 > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv;
  rm $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list_temp.tsv;

  # Look for insertions compared to database species
  cat $WORKDIR1/D4_align_shared_genes/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta;
  cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
  head -n2 | tail -n1 | grep -ob "[-\!]" | cut -f1 | awk '{print $0+1}' \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt;

  # Seperate gene into CDS sequences
  POS_START=1;
  for EXON_NR in $(seq 1 1 $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | wc -l)); do
    DB_SCAFFOLD=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f1);
    DB_START=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f2);
    DB_END=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f3);
    DB_STRANDEDNESS=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f4);
    MISSING1=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f9);
    MISSING2=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f10)
    EXON=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_exon_list.tsv | head -n $EXON_NR | tail -n1 | cut -f17);

    # Find end of CDS sequence
    POS_END=$(echo $POS_START + $DB_END - $DB_START | bc);

    # Look for insertions relative to database species and extend if necessary to end position
    X=$POS_START;
    Y=$POS_END;
    Z=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt | \
    awk -v X=$X -v Y=$Y 'BEGIN{sum=0} {if($1>=X && $1<=Y) {sum+=1}} END {print sum}');
    while [ $Z -gt 0 ]; do
      X=$(echo $Y + 1 | bc);
      Y=$(echo $Y + $Z | bc);
      Z=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_insertions.txt | \
      awk -v X=$X -v Y=$Y 'BEGIN{sum=0} {if($1>=X && $1<=Y) {sum+=1}} END {print sum}');
    done;
    POS_END=$Y;

    if [ $MISSING1 == "callable" ] && [ $MISSING2 == "callable" ]; then

      mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON;
      cd $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON;

      # Extract CDS sequence
      if [ $DB_STRANDEDNESS == "+" ]; then

        # Extract CDS sequence and add dummy bases in start of alignment
        seqkit subseq -r $POS_START:$POS_END \
        $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
        awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | \
        awk '{if($0 ~ /^>/) print; else print "XXXXXXXXXXXXXXXXXXXX"$0}' | tail -n+2 | tr '\!' '-' \
        > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_all_leadingdummy.fasta;

      else

        # Extract CDS sequence, add dummy bases in start of alignment
        seqkit subseq -r $POS_START:$POS_END \
        $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta | \
        awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | \
        while read alignment; do
          if [[ $alignment =~ ^'>' ]]; then
            echo $alignment;
          else
            echo $alignment | tr ACGTacgt TGCAtgca | rev;
          fi;
        done | awk '{if($0 ~ /^>/) print; else print "XXXXXXXXXXXXXXXXXXXX"$0}' | tail -n+2 | tr '\!' '-' \
        > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_all_leadingdummy.fasta;

      fi;

      # Convert alignment to vcf
      python $msa2vcf $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_all_leadingdummy.fasta | \
      awk -F'\t' -v DB_SCAFFOLD=$DB_SCAFFOLD -v DB_START=$DB_START \
      '{OFS="\t"; if($0 ~ /^##contig/) print "##contig=<ID="DB_SCAFFOLD">"; else if($0 ~ /^#/) print; \
      else if($2<=20) next; else {$1=DB_SCAFFOLD; $2=$2+DB_START-20-1; print}}' | \
      bcftools view -S $OUTDIR/metadata/individuals_$SUFFIX\_sampleID.txt | \
      bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON.vcf.gz;

      # Remove temporary files
      rm $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/$PROJECT1\_$PROJECT2\_$EXON\_alignment_NT_snpeff_all_leadingdummy.fasta \
      $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$EXON/snptable.txt;

      cd -;

   fi;

    # Update start position for next CDS sequence
    POS_START=$(echo $POS_END + 1 | bc);

  done;

  # Merge CDS back to gene
  find $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE -name "*.vcf.gz" \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_vcf_list.txt;
  bcftools concat -f $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_vcf_list.txt -Ou | \
  bcftools sort -T $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/temp_merge | \
  vcfallelicprimitives --keep-info --keep-geno | \
  awk -F $'\t' 'BEGIN {OFS = FS} /^[#]/ {print; next} {for (i=10; i<=NF; i++) { gsub("\\|","/",$i)} print}' | \
  bcftools norm -m -any --multi-overlaps "." | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  awk -F'\t' '{if($0 ~ /^#/) {print} else if($0 !~ /^#/ && $5 !~ /N/) {print}}' | \
  bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz;

  # Find alternate alleles present in outgroup lineages
  bcftools view $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz \
  -s$(echo $OUTGROUPS | tr '|' ',') | \
  bcftools norm -m +any | \
  vcffixup - | \
  vcffilter -f "AC > 0 " | \
  bgzip -c > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_outgroups.vcf.gz;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_outgroups.vcf.gz;

  ## Convert genotypes to haploid
  python $diploid2haploid \
  -i $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp1.vcf.gz \
  -o $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf;
  bgzip $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf;
  tabix $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE.vcf.gz;

  # Remove temporary files
  rm -r $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$GENE\_* \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_snpeff_all.fasta* \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_temp*;

};

## Excecute function in parallell
export -f shared_to_DB_ref;
parallel --colsep '\t' 'shared_to_DB_ref {}' :::: $GENES;

## Merge all genes across regions and find substitutions fixed across´lineages
for REGION1 in $(echo autosomal sex_phase sex_dropout); do

  if [ $REGION1 == "autosomal" ]; then

    REGION2=autosomal;
    find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" | grep -v "_outgroups.vcf.gz" \
    > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
    bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
    bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
    bcftools norm -m +any | \
    vcffixup - | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
    tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;

    bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz \
    -s^$(echo $DB_NAME"|"$OUTGROUPS | tr '|' ',') | \
    vcffixup - | \
    vcffilter -f "( AC / NS ) = 1" | \
    bcftools view -H | cut -f1,2 \
    > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.txt;
    bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz \
    -T $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.txt -s $DB_NAME | \
    bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;
    tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;

    bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
    > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
    >> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
    tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

    # Create lineage files
    for PROJECT in $(echo $PROJECT1 $PROJECT2); do

        if [ $PROJECT == $PROJECT1 ]; then
          REGION=$REGION1;
        else
          REGION=$REGION2;
        fi;

      bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
      -S $OUTDIR/metadata/$PROJECT\_individuals_A_sampleID.txt | \
      vcfallelicprimitives --keep-info --keep-geno | \
      vcffixup - | \
      vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_A.vcf.gz;
      tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_A.vcf.gz;

    done;

  else

    for REGION2 in $(echo sex_phase sex_dropout); do

      if [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_phase" ]; then

        find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" | grep -v "_outgroups.vcf.gz" \
        > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
        bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
        bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
        bcftools norm -m +any | \
        vcffixup - | \
        vcffilter -f "AC > 0" | \
        bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
        tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
  
        bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz \
        -s^$(echo $DB_NAME"|"$OUTGROUPS | tr '|' ',') | \
        vcffixup - | \
        vcffilter -f "( AC / NS ) = 1" | \
        bcftools view -H | cut -f1,2 \
        > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.txt;
        bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz \
        -T $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.txt -s $DB_NAME | \
        bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;
        tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;

      else

        find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" | grep -v "_outgroups.vcf.gz" \
        > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
        bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
        bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
        bcftools norm -m +any | \
        vcffixup - | \
        vcffilter -f "AC > 0" | \
        bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
        tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;

      fi;

      bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
      > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
      >> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
      tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

      # Create lineage files
      for PROJECT in $(echo $PROJECT1 $PROJECT2); do

        if [ $PROJECT == $PROJECT1 ]; then
            REGION=$REGION1;
        else
            REGION=$REGION2;
        fi;

        bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
        -S $OUTDIR/metadata/$PROJECT\_individuals_Z_sampleID.txt | \
        vcfallelicprimitives --keep-info --keep-geno | \
        vcffixup - | \
        vcffilter -f "AC > 0" | \
        bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_Z.vcf.gz;
        tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_Z.vcf.gz;

        if [ $REGION == "sex_phase" ]; then

          bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
          -S $OUTDIR/metadata/$PROJECT\_individuals_W_sampleID.txt | \
          vcfallelicprimitives --keep-info --keep-geno | \
          vcffixup - | \
          vcffilter -f "AC > 0" | \
          bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_W.vcf.gz;
          tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT\_$REGION\_W.vcf.gz;

        fi;

      done;

    done;

  fi;

done;

REGION1=sex_phase;
REGION2=autosomal;

find $OUTDIR/shared_genes/$REGION1\_$REGION2 -name "*.vcf.gz" | grep -v "_outgroups.vcf.gz" \
> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt;
bcftools concat -a -f $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_vcf_list.txt | \
bcftools sort --temp-dir $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp | \
bcftools norm -m +any | \
vcffixup - | \
vcffilter -f "AC > 0" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz;
bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz \
vcffilter -f "( AC / NS ) = 1" | \
bcftools view -s $DB_NAME | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_substitutions.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -h \
> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
bedtools sort -i $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp1.vcf.gz -faidx $DB_REF.fai \
>> $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
bgzip $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT1\_individuals_Z_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "AC > 0" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_Z.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_Z.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT1\_individuals_W_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "AC > 0" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_W.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT1\_$REGION1\_W.vcf.gz;

bcftools view $OUTDIR/merge_files/$REGION1\_$REGION2/$REGION1\_$REGION2\_temp2.vcf.gz \
-S $OUTDIR/metadata/$PROJECT2\_individuals_A_sampleID.txt | \
vcfallelicprimitives --keep-info --keep-geno | \
vcffixup - | \
vcffilter -f "AC > 0" | \
bgzip -c > $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT2\_$REGION2\_A.vcf.gz;
tabix $OUTDIR/merge_files/$REGION1\_$REGION2/$PROJECT2\_$REGION2\_A.vcf.gz;


# Merge lineage specific files
for REGION in $(echo A Z W); do

  for PROJECT in $(echo $PROJECT1 $PROJECT2); do

    find $OUTDIR/merge_files/ -name "$PROJECT\_*_$REGION.vcf.gz" \
    > $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt;

    bcftools concat -a -f $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt | \
    bcftools sort --temp-dir $OUTDIR/final_input/$PROJECT\_$REGION\_temp | \
    bcftools norm -m +any | \
    vcfallelicprimitives --keep-info --keep-geno | \
    vcffixup - | \
    vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz;
    tabix $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz;

    bcftools view $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz -h \
    > $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    bedtools sort -i $OUTDIR/final_input/$PROJECT\_$REGION\_temp1.vcf.gz -faidx $DB_REF.fai \
    >> $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    bgzip $OUTDIR/final_input/$PROJECT\_$REGION.vcf;
    tabix $OUTDIR/final_input/$PROJECT\_$REGION.vcf.gz;

    # Remove temporary files
    rm $OUTDIR/final_input/$PROJECT\_$REGION\_vcf_list.txt \
    $OUTDIR/final_input/$PROJECT\_$REGION\_temp*;

  done;

done;

# Merge outgroups
mkdir $OUTDIR/final_input/outgroups;
find $OUTDIR/shared_genes/ -name "*_outgroups.vcf.gz" \
> $OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups_vcf_list.txt;
bcftools concat -a -f $OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups_vcf_list.txt | \
bcftools sort --temp-dir $OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups_temp | \
bcftools norm -m +any | \
bgzip -c > $OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups.vcf.gz;
tabix $OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups.vcf.gz;

# Merge substitutions
mkdir $OUTDIR/final_input/substitutions;
find $OUTDIR/merge_files/ -name "*_substitutions.vcf.gz" \
> $OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions_vcf_list.txt;
bcftools concat -a -f $OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions_vcf_list.txt | \
bcftools sort --temp-dir $OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions_temp | \
bcftools norm -m +any | \
bgzip -c > $OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions.vcf.gz;
tabix $OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions.vcf.gz;

# Remove temporary files
rm -r $OUTDIR/shared_genes $OUTDIR/merge_files;

### Create database
cp $snpeff_config $OUTDIR/metadata/snpEff.config;
echo $DB_NAME genome, version $DB_NAME >> $OUTDIR/metadata/snpEff.config;
echo $DB_NAME.genome : $DB_NAME >> $OUTDIR/metadata/snpEff.config;

# Create transcripts only gtf
cat $GENES | cut -f2 > $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt;
cgat gtf2gtf --method=filter --filter-method=transcript --map-tsv-file $OUTDIR/metadata/$PROJECT1\_$PROJECT2\_transcripts.txt -I $DB_GTF \
-S $OUTDIR/data/$DB_NAME/genes.gtf;

# Copy files to the correct folders and names
cp $DB_REF $OUTDIR/data/genomes/$DB_NAME.fa;
cd $OUTDIR;

# Build database
snpEff build -gtf22 -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -treatAllAsProteinCoding -v $DB_NAME;


### Annotatate individual alleles
snpEff_inds() {

  IND=$1;
  SEX=$3;

  mkdir $OUTDIR/$PROJECT\_$REGION/$IND;
  cd $OUTDIR/$PROJECT\_$REGION/$IND;

  if [ $REGION == "autosomal" ]; then
    REGION1="A"
    REGION2="A"
    HAP1=A1;
    HAP2=A2;
  elif [ $REGION == "sexchrom" ]; then
    if [ $SEX == "Female" ]; then
      REGION1="Z"
      REGION2="W"
      HAP1=Z;
      HAP2=W;
    elif [ $SEX == "Male" ]; then
      REGION1="Z"
      REGION2="Z"
      HAP1=Z1;
      HAP2=Z2;
    fi;
  fi;

  # Subset to individual haplotype and remove across-lineage substitutions
  bcftools view $OUTDIR/final_input/$PROJECT\_$REGION1.vcf.gz -s $IND\_$HAP1 \
  -T ^$OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions.vcf.gz | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_nosub.vcf.gz;
  bcftools view $OUTDIR/final_input/$PROJECT\_$REGION2.vcf.gz -s $IND\_$HAP2 \
  -T ^$OUTDIR/final_input/substitutions/$PROJECT1\_$PROJECT2\_substitutions.vcf.gz | \
  vcffixup - | \
  vcffilter -f "AC > 0" | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_nosub.vcf.gz;

  # Estimate impact per allele for functional heterozygosity estimation
  for HAP in $(echo $HAP1 $HAP2); do

    snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_hets.html \
    -csvStats $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_hets.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_nosub.vcf.gz \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_hets1.vcf;
    bgzip $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_hets1.vcf;
    tabix $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_hets1.vcf.gz;

  done;

  # Remove variation segregating in outgroups
  bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_nosub.vcf.gz \
  -T ^$OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups.vcf.gz | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_nosub_noout.vcf.gz;
  bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_nosub.vcf.gz \
  -T ^$OUTDIR/final_input/outgroups/$PROJECT1\_$PROJECT2\_outgroups.vcf.gz | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_nosub_noout.vcf.gz;

  # Estimate impact per allele
  for HAP in $(echo $HAP1 $HAP2); do

    snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP.html \
    -csvStats $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_nosub_noout.vcf.gz \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_impacts.vcf;
    bgzip $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_impacts.vcf;
    tabix $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_impacts.vcf.gz;

    # Extract information from impact annotation
    bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_impacts.vcf.gz -H | \
    cut -f8 | cut -d'|' -f4,16 | cut -d',' -f1 | cut -d';' -f1 | grep WARNING | tr '|' '\t' | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_warnings.tsv;
    cat $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP.genes.txt | \
    awk -F'\t' -v PROJECT=$PROJECT -v HAP=$HAP 'BEGIN {COL=0} {if(NR==2) {for(i=1; i<=NF; i++) {if($i=="variants_impact_HIGH") COL=i}} \
    else if(COL > 0 && $COL > 0) {print $1"\t"$3"\t"$COL}}' \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_high_impacts.tsv;
    bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_impacts.vcf.gz -H | \
    cut -f8 | grep LOF | awk -F'LOF=\\(' '{if(NF==2) {print $2} else {for (i=2; i<NF; i++) printf $i "\n"; print $NF""}}' | \
    awk -F ')' '{print $1}' | cut -d'|' -f1 | tr -s '\n' | sort -u | awk -F'\t' -v PROJECT1=$PROJECT -v $HAP=HAP '{print $1"\t"1}' \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP\_LOF.tsv;

  done;

  # Estimate functional heterozygosity
  if [ $SEX == "Female" ]; then
    cat $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_LOF.tsv \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_LOF.tsv \
    $OUTDIR/metadata/$PROJECT\_dropout.txt | cut -f1 | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_LOF.txt;
    cat $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_high_impacts.tsv \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_high_impacts.tsv \
    $OUTDIR/metadata/$PROJECT\_dropout.txt | cut -f1 | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_high_impacts.txt;

  elif [ $SEX == "Male" ]; then
    cat $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_LOF.tsv \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_LOF.tsv | cut -f1 | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_LOF.txt;
    cat $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_high_impacts.tsv \
    $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_high_impacts.tsv | cut -f1 | sort -u \
    > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_high_impacts.txt;
  fi;

  bcftools merge $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_hets1.vcf.gz \
  $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_hets1.vcf.gz -m all --missing-to-ref | \
  grep -vFf $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_LOF.txt | \
  awk '{if($0 ~ /^#/ || $0 ~ /MODERATE/) print}' | \
  vcffixup - | \
  vcffilter -f "AC = 1" | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_LOF.vcf.gz;
  tabix $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_LOF.vcf.gz;

  bcftools merge $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP1\_hets1.vcf.gz \
  $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_$HAP2\_hets1.vcf.gz -m all --missing-to-ref | \
  grep -vFf $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_non_functional_high_impacts.txt | \
  awk '{if($0 ~ /^#/ || $0 ~ /MODERATE/) print}' | \
  vcffixup - | \
  vcffilter -f "AC = 1" | \
  bgzip -c > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_high.vcf.gz;
  tabix $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_high.vcf.gz;
  
  bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_LOF.vcf.gz -H | \
  awk -F'[\t;|]' '{print $14"\t"$13}' \
  > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_HETS_LOF.tsv;
    bcftools view $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_hets_MODERATE_high.vcf.gz -H | \
  awk -F'[\t;|]' '{print $14"\t"$13}' \
  > $OUTDIR/$PROJECT\_$REGION/$IND/$IND\_HETS_high.tsv;

  # Remove temporary files
  rm $(find $OUTDIR/$PROJECT\_$REGION/$IND -name "*.vcf.gz*" | grep -v "_impacts.vcf.gz" | grep -v MODERATE);

};

export -f snpEff_inds;

for PROJECT in $(echo $PROJECT1 $PROJECT2); do

  for REGION in $(echo autosomal sexchrom); do

    export PROJECT REGION;

    ## Excecute function in parallell
    parallel --header : --colsep '\t' 'snpEff_inds {}' :::: $MAINDIR/data/$PROJECT/metadata/sample_info.txt;

  done;

done;

export MAINDIR;

### Annotatate individual alleles
ind_geno() {

  GENE=$1;
  TRANS=$2;

  if [ $PROJECT == $PROJECT1 ]; then
    REGION=$3;
    STRATA=$5;
  elif [ $PROJECT == $PROJECT2 ]; then
    REGION=$4;
    STRATA=$6;
  fi;

  if [ $REGION != "autosomal" ]; then
    REGION0="sexchrom";
  else
    REGION0=$REGION;
  fi;

  for i in $(seq 1 1 $(cat $MAINDIR/data/$PROJECT/metadata/sample_info.txt | tail -n+2 | wc -l)); do
    IND=$(cat $MAINDIR/data/$PROJECT/metadata/sample_info.txt | tail -n+2 | head -n$i | tail -n1 | cut -f1);
    SEX=$(cat $MAINDIR/data/$PROJECT/metadata/sample_info.txt | tail -n+2 | head -n$i | tail -n1 | cut -f3);

    if [ $REGION == "autosomal" ]; then
      REGION1="A";
      REGION2="A";
      HAP1=A1;
      HAP2=A2;
      DROPOUT2=0;
      LOF1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
      HIGH1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');
      LOF2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
      HIGH2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');

    elif [ $REGION == "sex_phase" ] || [ $REGION == "partial_dropout" ] || [ $REGION == "sex_dropout" ]; then

      if [ $SEX == "Female" ]; then
        REGION1="Z"
        REGION2="W"
        HAP1=Z;
        HAP2=W;
        LOF1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
        HIGH1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');
        if [ $REGION == "sex_dropout" ] || [ $REGION == "partial_dropout" ]; then
          DROPOUT2=1;
          LOF2=0;
          HIGH2=0;
        elif [ $REGION == "sex_phase" ]; then
          DROPOUT2=0;
          LOF2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
          HIGH2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');
        fi;

      elif [ $SEX == "Male" ]; then
        REGION1="Z"
        REGION2="Z"
        HAP1=Z1;
        HAP2=Z2;
        DROPOUT2=0;
        LOF1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
        HIGH1=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP1\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');
        LOF2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
        HIGH2=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_$HAP2\_high_impacts.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=$3; print SUM; exit}} END {if(SUM==0) print 0}');

      fi;

    fi;

    DOM_SCORE_LOF=$(echo $LOF1 $LOF2 $DROPOUT2 | awk '{if($1==0 || ($2==0 && $3==0)) {SUM=0} else {SUM=1}; print SUM}');
    ADD_SCORE_LOF=$(echo $LOF1 $LOF2 $DROPOUT2 | awk '{if($2==1 || $3==1) {HAP2=1} else (HAP2=0); print ($1+HAP2)/2}');
    if (( $(echo "$ADD_SCORE_LOF == 0" | bc -l) )); then
      HET_SCORE_LOF=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_HETS_LOF.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
      SHELTER_HAP_LOF=NA;
    else
      HET_SCORE_LOF=0;
      if (( $(echo "$ADD_SCORE_LOF < 1" | bc -l) )) && [ $REGION != "autosomal" ] && [ $SEX == "Female" ]; then
        SHELTER_HAP_LOF=$(echo $LOF1 $LOF2 $DROPOUT2 | awk '{if($2==1 || $3==1) {HAP2=1} else (HAP2=0); if($1 > HAP2) {HAP="W"} else {HAP="Z"}; print HAP}');
      else
        SHELTER_HAP_LOF=NA;
      fi;
    fi;
    DOM_SCORE_HIGH=$(echo $HIGH1 $HIGH2 $DROPOUT2 | awk '{if($1==0 || ($2==0 && $3==0)) {SUM=0} else {SUM=1}; print SUM}');
    ADD_SCORE_HIGH=$(echo $HIGH1 $HIGH2 $DROPOUT2 | awk '{if($1 > 0) {$1=1}; if($2 > 0 || $3 > 0) {HAP2=1} else {HAP2=0}; print ($1+HAP2)/2}');
    if (( $(echo "$ADD_SCORE_HIGH == 0" | bc -l) )); then
      HET_SCORE_HIGH=$(cat $OUTDIR/$PROJECT\_$REGION0/$IND/$IND\_HETS_high.tsv | awk -F'\t' -v GENE=$GENE 'BEGIN {SUM=0} {if($1==GENE) {SUM=1; print SUM; exit}} END {if(SUM==0) print 0}');
      SHELTER_HAP_HIGH=NA;
    else
      HET_SCORE_HIGH=0;
      if (( $(echo "$ADD_SCORE_HIGH < 1" | bc -l) )) && [ $REGION != "autosomal" ] && [ $SEX == "Female" ]; then
        SHELTER_HAP_HIGH=$(echo $HIGH1 $HIGH2 $DROPOUT2 | awk '{if($2 > 0 || $3 > 0) {HAP2=1} else (HAP2=0); if($1 > HAP2) {HAP="W"} else {HAP="Z"}; print HAP}');
      else
        SHELTER_HAP_HIGH=NA;
      fi;
    fi;

    # Report data
    echo -e "${GENE}\t${TRANS}\t${REGION}\t${STRATA}\t${PROJECT}\t${IND}\t${SEX}\t${REGION1}\t${LOF1}\t${HIGH1}\t${REGION2}\t${DROPOUT2}\t${LOF2}\t${HIGH2}\t${DOM_SCORE_LOF}\t${ADD_SCORE_LOF}\t${HET_SCORE_LOF}\t${DOM_SCORE_HIGH}\t${ADD_SCORE_HIGH}\t${HET_SCORE_HIGH}\t${SHELTER_HAP_LOF}\t${SHELTER_HAP_HIGH}";

  done > $OUTDIR/genotype_info/$PROJECT\_$GENE\_info.tsv;

};

export -f ind_geno;

for PROJECT in $(echo $PROJECT1 $PROJECT2); do

  export PROJECT;

  ## Excecute function in parallell
  parallel --header : --colsep '\t' 'ind_geno {}' :::: $GENES;

done;

# Merge all gene info
echo "geneID transcriptID Region Strata Project Individual Sex Hap1_region LOF1 HIGH1 Hap2_region W_Dropout LOF2 High2 Dominance_score_LOF Additive_score_LOF Heterozygosity_score_LOF Dominance_score_HIGH Additive_score_HIGH Heterozygosity_score_HIGH Sheltering_gametolog_LOF Sheltering_gametolog_High" | \
tr ' ' '\t' \
> $OUTDIR/metadata/genotype_info.txt;
cat $(find $OUTDIR/genotype_info -name "*_info.tsv") | sort -k1,1 -k5,5 -k7,7 -k6,6 \
>> $OUTDIR/metadata/genotype_info.txt;

# Remove temporary files
rm -r $OUTDIR/genotype_info;

# Annotate variants with outgroups and all substitutions
filter_modifier() {
  bcftools view $1 | \
  awk 'BEGIN {OFS="\t"} {if ($0 ~ /^#/ ) {print $0} else { 
      split($8, info, ";"); 
      for (i in info) { 
          if (info[i] ~ /ANN=/) { 
              ann = substr(info[i], 5);  # Skip the "ANN=" prefix
              split(ann, effects, ",");  # Split the annotations into individual effects
              new_ann = ""; 
              for (j in effects) { 
                  if (effects[j] !~ /MODIFIER/) { 
                      if (new_ann != "") new_ann = new_ann "," effects[j]; 
                      else new_ann = effects[j];}}
              if (new_ann != "") {
                  info[i] = "ANN=" new_ann;
              } else {
                  info[i] = "";}}}
      $8 = info[1];
      for (i = 2; i <= length(info); i++) {
          if (info[i] != "") {
              $8 = $8 ";" info[i];}}
      print $0;}}' | \
  bgzip -c > $(echo $1 | rev | cut -d'.' -f3- | rev | awk '{print $0"_no_MODIFIER.vcf.gz"}');
  tabix $(echo $1 | rev | cut -d'.' -f3- | rev | awk '{print $0"_no_MODIFIER.vcf.gz"}');
};

bcftools merge $OUTDIR/final_input/$PROJECT1\_A.vcf.gz \
$OUTDIR/final_input/$PROJECT2\_A.vcf.gz --missing-to-ref --force-samples -m all | \
bcftools view -s ^$(echo $DB_NAME"|"$OUTGROUPS | tr '|' '\n' | awk '{print "2:"$1","}' | tr -d '\n' | sed 's/.$//') | \
bgzip > $OUTDIR/final_input/$PROJECT1\_$PROJECT2\_A.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A_impacts.html \
-csvStats $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT1\_$PROJECT2\_A.vcf.gz \
> $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A_impacts.vcf;
bgzip $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A_impacts.vcf;
tabix $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A_impacts.vcf.gz;
filter_modifier $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_A_impacts.vcf.gz;

bcftools merge $OUTDIR/final_input/$PROJECT1\_Z.vcf.gz \
$OUTDIR/final_input/$PROJECT2\_Z.vcf.gz --missing-to-ref --force-samples -m all | \
bcftools view -s ^$(echo $DB_NAME"|"$OUTGROUPS | tr '|' '\n' | awk '{print "2:"$1","}' | tr -d '\n' | sed 's/.$//') | \
bgzip > $OUTDIR/final_input/$PROJECT1\_$PROJECT2\_Z.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z_impacts.html \
-csvStats $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT1\_$PROJECT2\_Z.vcf.gz \
> $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z_impacts.vcf;
bgzip $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z_impacts.vcf;
tabix $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z_impacts.vcf.gz;
filter_modifier $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_Z_impacts.vcf.gz;

bcftools merge $OUTDIR/final_input/$PROJECT1\_W.vcf.gz \
$OUTDIR/final_input/$PROJECT2\_W.vcf.gz --missing-to-ref --force-samples -m all | \
bcftools view -s ^$(echo $DB_NAME"|"$OUTGROUPS | tr '|' '\n' | awk '{print "2:"$1","}' | tr -d '\n' | sed 's/.$//') | \
bgzip > $OUTDIR/final_input/$PROJECT1\_$PROJECT2\_W.vcf.gz;

snpEff -c $OUTDIR/metadata/snpEff.config -dataDir $OUTDIR/data -s $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W_impacts.html \
-csvStats $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W.csv -treatAllAsProteinCoding -v -d -lof $DB_NAME \
$OUTDIR/final_input/$PROJECT1\_$PROJECT2\_W.vcf.gz \
> $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W_impacts.vcf;
bgzip $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W_impacts.vcf;
tabix $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W_impacts.vcf.gz;
filter_modifier $OUTDIR/all_data/$PROJECT1\_$PROJECT2\_W_impacts.vcf.gz;
