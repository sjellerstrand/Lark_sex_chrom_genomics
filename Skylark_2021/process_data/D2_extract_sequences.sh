#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 50:00:00
#SBATCH -J extract_sequences
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
FEATURES=$WORKDIR/D1_filter_genes/$PROJECT\_features.tsv;
VCFS=$WORKDIR/A8_PhaseWY/final_output/vcfs;
BEDS=$WORKDIR/A8_PhaseWY/final_output/beds;
REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
DB_NAME=Tgut;
DB_REF=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.fasta;
METADATA=$MAINDIR/data/$PROJECT/metadata;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");

### Load modules
module load bioinfo-tools samtools/1.14 bcftools/1.14 vcflib/1.0.1 SeqKit/0.15.0 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/D2_extract_sequences;
OUTDIR=$OUTDIR/D2_extract_sequences;
mkdir $OUTDIR/CDS \
$OUTDIR/genes \
$OUTDIR/genes/autosomal \
$OUTDIR/genes/sex_dropout \
$OUTDIR/genes/partial_dropout \
$OUTDIR/genes/sex_phase \
$OUTDIR/metadata;

# Create outgroup list
echo $OUTGROUPS | tr ":" "\n" | tr " " "\t" \
> $OUTDIR/metadata/outgroups.txt;

export OUTDIR FEATURES PROJECT METADATA VCFS BEDS REF DB_NAME DB_REF;

extract_sequences() {
  DB_SCAFFOLD=$1;
  DB_START=$2;
  DB_END=$3;
  DB_STRANDEDNESS=$4;
  REF_SCAFFOLD=$5;
  REF_START=$6;
  REF_END=$7;
  REF_STRANDEDNESS=$8;
  MISSING=$9;
  SHARED=${10};
  REGION=${11};
  STRATA=${12};
  EXON=${13};

  echo $EXON
  mkdir $OUTDIR/CDS/$EXON;

  # Data type specific settings
  if [ $REGION == "autosomal" ]; then
    FEMALE_SUFFIX1=$(echo A1);
    FEMALE_SUFFIX2=$(echo A2);
    MALE_SUFFIX1=$(echo A1);
    MALE_SUFFIX2=$(echo A2);
    VCF_REGION1=$(echo autosomal);
    VCF_REGION2=$(echo autosomal);
    PHASE1=$(echo 1pIu);
    PHASE2=$(echo 2pIu);

  elif [ $REGION == "sex_dropout" ]; then
    FEMALE_SUFFIX2=$(echo Z);
    MALE_SUFFIX1=$(echo Z1);
    MALE_SUFFIX2=$(echo Z2);
    VCF_REGION2=$(echo homogametic);
    PHASE2=$(echo 1);

  elif [ $REGION == "sex_phase" ] || [ $REGION == "partial_dropout" ] ; then
    FEMALE_SUFFIX1=$(echo W);
    FEMALE_SUFFIX2=$(echo Z);
    MALE_SUFFIX1=$(echo Z1);
    MALE_SUFFIX2=$(echo Z2);
    VCF_REGION1=$(echo heterogametic);
    VCF_REGION2=$(echo homogametic);
    PHASE1=$(echo 1);
    PHASE2=$(echo 1);
  fi;

  # Extract region from database reference
  echo $DB_SCAFFOLD $DB_START $DB_END | \
  awk '{print $1":"$2"-"$3}' \
  > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_DB_ref.txt;
  echo -e ">${DB_NAME}" \
  > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

  if [ $DB_STRANDEDNESS == "+" ]; then
    samtools faidx $DB_REF -r $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_DB_ref.txt | \
    grep -v "^>" | tr -d "\n" | tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
    >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
  else
    samtools faidx $DB_REF -r $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_DB_ref.txt | \
    grep -v "^>" | tr -d "\n" | tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
    tr ACGTacgt TGCAtgca | rev \
    >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
  fi;

  # Extract region from outgroup reference
  for OUT in $(cat $OUTDIR/metadata/outgroups.txt | cut -f1); do
    echo -e ">${OUT}" \
    >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

    if [ $MISSING != "callable" ] || [ $SHARED != "shared" ]; then
      DB_LENGTH=$(echo $DB_END - $DB_START + 1 | bc);
      awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

    else

      EXON_INFO=$(cat $(cat $OUTDIR/metadata/outgroups.txt | grep $OUT | cut -f 3) | grep $EXON | cut -f1,4,5,7);
      if [ $(echo $EXON_INFO | tr -d '\n' | wc -c) -gt 0 ]; then

        echo $EXON_INFO | awk '{print $1":"$2"-"$3}' \
        > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_$OUT.txt;

        if [ $(echo $EXON_INFO | cut -d' ' -f4) == "+" ]; then
          samtools faidx $(cat $OUTDIR/metadata/outgroups.txt | grep $OUT | cut -f2) \
          -r $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_$OUT.txt | \
          grep -v "^>" | tr -d "\n" | tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        else
          samtools faidx $(cat $OUTDIR/metadata/outgroups.txt | grep $OUT | cut -f2) -r \
          $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_$OUT.txt | \
          grep -v "^>" | tr -d "\n" | tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
          tr ACGTacgt TGCAtgca | rev \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        fi;

      else
        DB_LENGTH=$(echo $DB_END - $DB_START + 1 | bc)
        awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      fi;

    fi;

  done;

  if [ $MISSING != "callable" ]; then

    # Add females
    for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do

      # Add first or W-haplotype
      if [ $REGION != "sex_dropout" ]; then
        echo -e ">${IND}_${FEMALE_SUFFIX1}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      fi;

      # Add second or Z-haplotype
      echo -e ">${IND}_${FEMALE_SUFFIX2}" \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
    done;

    # Add males
    for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do

      # Add first haplotype
      echo -e ">${IND}_${MALE_SUFFIX1}" \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

      # Add second haplotype
      echo -e ">${IND}_${MALE_SUFFIX2}" \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      awk -v DB_LENGTH=$DB_LENGTH 'BEGIN {OFS="N"; $DB_LENGTH="N"; print}' \
      >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
    done;

    # Remove temporary files
    rm $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_*.txt;

  else

    # Extract region from species reference
    echo $REF_SCAFFOLD $REF_START $REF_END | \
    awk '{print $1":"$2"-"$3}' \
    > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.txt;
    samtools faidx $REF -r $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.txt \
    > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta;

    # Combine overlapping alleles
    bcftools view $VCFS/$PROJECT\_$VCF_REGION2.vcf.gz -r $(cat $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.txt) | \
    vcfcreatemulti | vcffixup - | vcffilter -f "AC > 0" | \
    bgzip -c > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_2.vcf.gz;
    tabix $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_2.vcf.gz;
    VCF2=$OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_2.vcf.gz;

    if [ $REGION != "sex_phase" ] && [ $REGION != "partial_dropout" ]; then
      VCF1=$OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_2.vcf.gz;
    else
      bcftools view $VCFS/$PROJECT\_$VCF_REGION1.vcf.gz -r $(cat $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.txt) | \
      vcfcreatemulti | vcffixup - | vcffilter -f "AC > 0" | \
      bgzip -c > $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_1.vcf.gz;
      tabix $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_1.vcf.gz;
      VCF1=$OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_1.vcf.gz;
    fi;

    if [ $REF_STRANDEDNESS == "+" ]; then

      # Add females
      for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do

        if [ $REGION != "autosomal" ]; then

          if [ $REGION == "sex_phase" ] || [ $REGION == "partial_dropout" ]; then


            # Remove deletions extending into start of sequence
            DIFF_START1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START '{if($1 < REF_START) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START1=1;
            fi;

            # Remove insertion extending from end of sequence
            ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
            awk -F'\t' '{print $NF}' | cut -d'|' -f1);
            if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_END1=1;
            fi;
            DIFF_END1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
            if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_END1=1;
            fi;

          fi;

          # Remove deletions extending into start of sequence
          DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START '{if($1 < REF_START) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          fi;

          # Remove insertion extending from end of sequence
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;
          DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;

        else

          # Remove deletions extending into start of sequence
          ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START1=1;
          else
            DIFF_START1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START1=1;
            fi;
          fi;
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f2);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          else
            DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START2=1;
            fi;
          fi;

          # Remove insertion extending from end of sequence
          ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END1=1;
          fi;
          DIFF_END1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END1=1;
          fi;
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f2);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;
          DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;

        fi;

        # Add first or W-haplotype
        if [ $REGION != "sex_dropout" ]; then
          echo -e ">${IND}_${FEMALE_SUFFIX1}" \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
          samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
          bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H $PHASE1 $VCF1 \
          -m $BEDS/$PROJECT\_missing_region.bed | \
          grep -v "^>" | tr -d "\n" | cut -c$DIFF_START1- | rev | cut -c$DIFF_END1- | rev | \
          tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        fi;

        # Add second or Z-haplotype
        echo -e ">${IND}_${FEMALE_SUFFIX2}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H $PHASE2 $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START2- | rev | cut -c$DIFF_END2- | rev | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      done;

      # Add males
      for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do

        # Remove deletions extending into start of sequence
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f1);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_START1=1;
        else
          DIFF_START1=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START1=1;
          fi;
        fi;
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f2);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_START2=1;
        else
          DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          fi;
        fi;

        # Remove insertion extending from end of sequence
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f1);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END1=1;
        fi;
        DIFF_END1=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
        awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
        if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END1=1;
        fi;
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f2);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END2=1;
        fi;
        DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
        awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
        if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END2=1;
        fi;

        # Add first haplotype
        echo -e ">${IND}_${MALE_SUFFIX1}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H 1pIu $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START1- | rev | cut -c$DIFF_END1- | rev | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

        # Add second haplotype
        echo -e ">${IND}_${MALE_SUFFIX2}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H 2pIu $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START2- | rev | cut -c$DIFF_END2- | rev | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      done;

    else

      # Add females
      for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Female") print $1}'); do

        if [ $REGION != "autosomal" ]; then

          if [ $REGION == "sex_phase" ] || [ $REGION == "partial_dropout" ]; then

            # Remove deletions extending into start of sequence
            DIFF_START1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START '{if($1 < REF_START) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START1=1;
            fi;

            # Remove insertion extending from end of sequence
            ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
            awk -F'\t' '{print $NF}' | cut -d'|' -f1);
            if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_END1=1;
            fi;
            DIFF_END1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
            if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_END1=1;
            fi;

          fi;

          # Remove deletions extending into start of sequence
          DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START '{if($1 < REF_START) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          fi;

          # Remove insertion extending from end of sequence
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;
          DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;

        else

          # Remove deletions extending into start of sequence
          ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START1=1;
          else
            DIFF_START1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START1=1;
            fi;
          fi;
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f2);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          else
            DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
            awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
            if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
              DIFF_START2=1;
            fi;
          fi;

          # Remove insertion extending from end of sequence
          ALLELE=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f1);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END1=1;
          fi;
          DIFF_END1=$(bcftools view $VCF1 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END1=1;
          fi;
          ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
          awk -F'\t' '{print $NF}' | cut -d'|' -f2);
          if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;
          DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
          if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_END2=1;
          fi;

        fi;

        # Add first or W-haplotype
        if [ $REGION != "sex_dropout" ]; then
          echo -e ">${IND}_${FEMALE_SUFFIX1}" \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
          samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
          bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H $PHASE1 $VCF1 \
          -m $BEDS/$PROJECT\_missing_region.bed | \
          grep -v "^>" | tr -d "\n" | cut -c$DIFF_START1- | \
          tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
          tr ACGTacgt TGCAtgca | rev | cut -c$DIFF_END1- \
          >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        fi;

        # Add second or Z-haplotype
        echo -e ">${IND}_${FEMALE_SUFFIX2}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H $PHASE2 $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START2- | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
        tr ACGTacgt TGCAtgca | rev | cut -c$DIFF_END2- \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      done;

      # Add males
      for IND in $(cat $METADATA/sample_info.txt | awk -F'\t' '{if($3=="Male") print $1}'); do

        # Remove deletions extending into start of sequence
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f1);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_START1=1;
        else
          DIFF_START1=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START1 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START1=1;
          fi;
        fi;
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f2);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_START2=1;
        else
          DIFF_START2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | head -n1 | cut -f2,4,5 | \
          awk -F'[\t,]' -v REF_START=$REF_START -v ALLELE=$ALLELE '{if($1 < REF_START && $2 != $(ALLELE+2)) {print $1+length($2)-REF_START+1}}');
          if [ $(echo $DIFF_START2 | tr -d '\n' | wc -c) -eq 0 ]; then
            DIFF_START2=1;
          fi;
        fi;

        # Remove insertion extending from end of sequence
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f1);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END1=1;
        fi;
        DIFF_END1=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
        awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
        if [ $(echo $DIFF_END1 | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END1=1;
        fi;
        ALLELE=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | \
        awk -F'\t' '{print $NF}' | cut -d'|' -f2);
        if [ $(echo $ALLELE | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END2=1;
        fi;
        DIFF_END2=$(bcftools view $VCF2 -s $IND | vcffilter -f "AC > 0" | bcftools view -H | tail -n1 | cut -f2,4,5 | \
        awk -F'[\t,]' -v REF_END=$REF_END -v ALLELE=$ALLELE '{if($1 == REF_END) {print $1+length($(ALLELE+2))-REF_END}}');
        if [ $(echo $DIFF_END2 | tr -d '\n' | wc -c) -eq 0 ]; then
          DIFF_END2=1;
        fi;

        # Add first haplotype
        echo -e ">${IND}_${MALE_SUFFIX1}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H 1pIu $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START1- | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
        tr ACGTacgt TGCAtgca | rev | cut -c$DIFF_END1- \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;

        # Add second haplotype
        echo -e ">${IND}_${MALE_SUFFIX2}" \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
        samtools faidx $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta | \
        bcftools consensus -s $IND -f $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta -H 2pIu $VCF2 \
        -m $BEDS/$PROJECT\_missing_region.bed | \
        grep -v "^>" | tr -d "\n" | cut -c$DIFF_START2- | \
        tr '[:lower:]' '[:upper:]' | awk '{print  $0"\n"}' | sed '/^$/d' | \
        tr ACGTacgt TGCAtgca | rev | cut -c$DIFF_END2- \
        >> $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta;
      done;

    fi;

    # Remove temporary files
    rm $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_*.txt \
    $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_region_species_ref.fasta*;

  fi;

};

## Excecute function in parallell
export -f extract_sequences;
parallel --colsep '\t' 'extract_sequences {}' :::: $FEATURES;

cat $FEATURES | cut -f11,12,13 | \
awk -F'[\t_]' '{if($1=="autosomal") {print $3"\t"$4"_"$5"\t"$1"\t"$2} \
else {print $4"\t"$5"_"$6"\t"$1"_"$2"\t"$3}}' | sort -u \
> $OUTDIR/metadata/genes.tsv;

concatenate_CDS() {

  GENE=$1;
  TRANS=$2;
  REGION=$3;
  STRATA=$4;

  echo $GENE;
  mkdir $OUTDIR/genes/$REGION/$GENE;

  # Sort CDS sequences by exon order
  for EXON in $(cat $FEATURES | cut -f13 | grep $TRANS); do
    echo $OUTDIR/CDS/$EXON/$PROJECT\_$EXON\_all_sequences.fasta $(echo $EXON | rev | cut -d'_' -f2 | rev);
  done | tr ' ' '\t' | sort -nk2,2 | cut -f1 > $OUTDIR/genes/$REGION/$GENE/$GENE\_exon_list.tsv;

  # Merge CDS sequences into full transcript sequence
  if [ $(cat $OUTDIR/genes/$REGION/$GENE/$GENE\_exon_list.tsv | wc -l) -gt 1 ]; then
    seqkit concat $(cat $OUTDIR/genes/$REGION/$GENE/$GENE\_exon_list.tsv) \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_all_sequences.fasta;
  else
    cp $(cat $OUTDIR/genes/$REGION/$GENE/$GENE\_exon_list.tsv) \
    $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_all_sequences.fasta;
  fi;

  # Split genes with W sequences into two seperate files
  if [ $REGION == "sex_phase" ] || [ $REGION == "partial_dropout" ]; then
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_all_sequences.fasta | \
    awk '{if($0 ~ ">" && $0 ~ /W/) W="T"; else if($0 ~ ">" && $0 !~ /W/) W="F"; if(W=="F") print}' \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_nonW_sequences.fasta;
    cat $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_all_sequences.fasta | \
    awk '{if($0 ~ ">" && $0 ~ /W/) W="T"; else if($0 ~ ">" && $0 !~ /W/) W="F"; if(W=="T") print}' \
    > $OUTDIR/genes/$REGION/$GENE/$PROJECT\_$GENE\_W_sequences.fasta;
  fi;

  # Remove temporary files
  rm -r $OUTDIR/CDS/$GENE\_$TRANS\_* \
  $OUTDIR/genes/$REGION/$GENE/$GENE\_exon_list.tsv;
};

## Excecute function in parallell
export -f concatenate_CDS;
parallel --colsep '\t' 'concatenate_CDS {}' :::: $OUTDIR/metadata/genes.tsv;

# Remove temporary files
rm -r $OUTDIR/CDS;
