#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J filter_genes
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Rasolark_2021;
PROJECT2=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
METADATA2=$MAINDIR/data/$PROJECT2/metadata;
BEDS1=$WORKDIR/A8_PhaseWY/final_output/beds;
BEDS2=$WORKDIR/A9_define_regions;
SPEC_NAME=Aarv;
DB_NAME=Tgut;
DB_REF=$WORKDIR2/B3_annotation_lift_over/$DB_NAME.fa;
SPEC_REF=$WORKDIR/A4_make_new_reference/$PROJECT\_consensus_reference.fasta;
SPEC_GTF=$WORKDIR2/B3_annotation_lift_over/$SPEC_NAME\_$DB_NAME\_liftover.gtf_polished;
DB_GTF=$WORKDIR2/B3_annotation_lift_over/$DB_NAME\_modified.gtf;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR2}/B3_annotation_lift_over_Pmaj/Pmaj_${DB_NAME}_liftover.gtf_polished
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR2}/B3_annotation_lift_over_Falb/Falb_${DB_NAME}_liftover.gtf_polished");

### Load modules
module load bioinfo-tools BEDTools/2.29.2;

### Activate conda environment
conda activate liftover;

### Create folders
mkdir $OUTDIR/D1_filter_genes;
OUTDIR=$OUTDIR/D1_filter_genes;

# Create unique identifiers for each feature
SPECIES=$(echo -e "${SPEC_NAME} ${SPEC_REF} ${SPEC_GTF}:${DB_NAME} ${DB_REF} ${DB_GTF}:${OUTGROUPS}");
for SPEC in $(echo $SPECIES | tr ":" "\n" | cut -d' ' -f1); do
  cat $(echo $SPECIES | tr ":" "\n" | awk -v SPEC=$SPEC '{if($1==SPEC) print $3}') | \
  cut -f3,9 | tr "\t" " " | \
  awk -F' ' '{ \
    if($0 ~ /#/) { \
      print \
    } \
    else if($1=="gene") { \
      for(i=2; i<=NF; i++) { \
        if($i ~ /^gene_id/) { \
          print $(i+1)"_parent_gene"; next
        } \
      } \
    }
    else if($1=="transcript") { \
      for(i=2; i<=NF; i++) { \
        if($i ~ /^gene_id/) { \
          gene=$(i+1) \
        } \
        if($i ~ /^transcript_id/) { \
          transcript=$(i+1) \
        } \
      } \
      print gene"_"transcript"_parent_transcript" \
    } \
    else { \
      for(i=2; i<=NF; i++) { \
        if($i ~ /^gene_id/) { \
          gene=$(i+1) \
        } \
        if($i ~ /^transcript_id/) { \
          transcript=$(i+1) \
        } \
        if($i ~ /^exon_number/) { \
          exon=$(i+1) \
        } \
      } \
      print gene"_"transcript"_exon_"exon"_"$1
    } \
  }' | \
  tr -d '"' | tr -d ';' \
  > $OUTDIR/$SPEC\_features_IDs.tsv;
  paste -d "\t" $(echo $SPECIES | tr ":" "\n" | awk -v SPEC=$SPEC '{if($1==SPEC) print $3}') $OUTDIR/$SPEC\_features_IDs.tsv \
  > $OUTDIR/$SPEC\_uniqID.gtf;
done;

## Remove non-unique feature IDs and corresponding transcripts
cat $OUTDIR/$SPEC_NAME\_uniqID.gtf | awk '{print $NF}' | sort | uniq -D | \
sed '/^$/d' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/duplicateIDs.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID.gtf | awk '{print $NF}' | grep -v "parent" | cut -d'_' -f1,2,3  | sort -u \
> $OUTDIR/transcripts_uniqfeats.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID.gtf | grep -vFf $OUTDIR/duplicateIDs.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats.gtf;

## Find transcripts in DB species with start and stop codon
cat $OUTDIR/$DB_NAME\_uniqID.gtf | \
awk -F'\t' '{if($3=="start_codon") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3  | sort -u \
> $OUTDIR/transcripts_$DB_NAME\_start_codon.txt;
cat $OUTDIR/$DB_NAME\_uniqID.gtf | \
awk -F'\t' '{if($3=="stop_codon") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3  | sort -u \
> $OUTDIR/transcripts_$DB_NAME\_stop_codon.txt;
cat $OUTDIR/transcripts_$DB_NAME\_start_codon.txt | grep -Ff $OUTDIR/transcripts_$DB_NAME\_stop_codon.txt \
> $OUTDIR/transcripts_$DB_NAME\_start_stop_codon.txt;
 cat $OUTDIR/transcripts_uniqfeats.txt | grep -vFf $OUTDIR/transcripts_$DB_NAME\_start_stop_codon.txt \
> $OUTDIR/transcripts_$DB_NAME\_no_start_stop_codon.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats.gtf | grep -vFf $OUTDIR/transcripts_$DB_NAME\_no_start_stop_codon.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop.gtf;

## Remove non-gene and non-transcript features not overlapping properly with callable regions
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop.gtf \
-b $BEDS2/$PROJECT\_target_region.bed -v -f 0.80  | \
awk -F'\t' '{if($0 !~ /^#/ && $3!="gene" && $3!="transcript") print}' | awk '{print $NF}' \
> $OUTDIR/uncallable_features.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop.gtf | grep -vFf $OUTDIR/uncallable_features.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf;

## Remove transcripts present in both autosomal and sex-linked regions
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/all_callable_transcripts.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf \
-b $BEDS1/$PROJECT\_autosomal.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_autosomal_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf \
-b $BEDS1/$PROJECT\_hetgam_dropout.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_hetgam_dropout_temp1.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf \
-b $BEDS1/$PROJECT\_homogametic.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_homogametic_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf \
-b $BEDS1/$PROJECT\_heterogametic.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_heterogametic_temp1.txt;

cat $OUTDIR/transcripts_autosomal_temp.txt | grep -vFf $OUTDIR/transcripts_homogametic_temp.txt \
> $OUTDIR/transcripts_autosomal.txt;
cat $OUTDIR/transcripts_homogametic_temp.txt | grep -vFf $OUTDIR/transcripts_autosomal_temp.txt \
> $OUTDIR/transcripts_homogametic.txt;
cat $OUTDIR/transcripts_homogametic.txt | grep -Ff $OUTDIR/transcripts_hetgam_dropout_temp1.txt \
> $OUTDIR/transcripts_hetgam_dropout_temp2.txt;
cat $OUTDIR/transcripts_homogametic.txt | grep -Ff $OUTDIR/transcripts_heterogametic_temp1.txt \
> $OUTDIR/transcripts_heterogametic_temp2.txt;
cat $OUTDIR/transcripts_hetgam_dropout_temp2.txt | grep -vFf $OUTDIR/transcripts_heterogametic_temp2.txt \
> $OUTDIR/transcripts_hetgam_dropout.txt;
cat $OUTDIR/transcripts_heterogametic_temp2.txt | grep -vFf $OUTDIR/transcripts_hetgam_dropout_temp2.txt \
> $OUTDIR/transcripts_heterogametic.txt;
cat $OUTDIR/transcripts_homogametic.txt | grep -vFf $OUTDIR/transcripts_hetgam_dropout.txt | \
grep -vFf $OUTDIR/transcripts_heterogametic.txt \
> $OUTDIR/transcripts_partial_hetgam_dropout.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf | \
grep -Ff $OUTDIR/transcripts_partial_hetgam_dropout.txt | \
bedtools intersect -a - -b $BEDS1/$PROJECT\_hetgam_dropout.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | sort -u \
> $OUTDIR/transcripts_partial_dropout_feature_dropout.txt;

cat $OUTDIR/transcripts_autosomal.txt | awk '{print $0"\tautosomal"}' \
> $OUTDIR/transcripts_all_nonconfunding_region1.txt;
cat $OUTDIR/transcripts_hetgam_dropout.txt | awk '{print $0"\tsex_dropout"}' \
>> $OUTDIR/transcripts_all_nonconfunding_region1.txt;
cat $OUTDIR/transcripts_heterogametic.txt | awk '{print $0"\tsex_phase"}' \
>> $OUTDIR/transcripts_all_nonconfunding_region1.txt;
cat $OUTDIR/transcripts_partial_hetgam_dropout.txt | awk '{print $0"\tpartial_dropout"}' \
>> $OUTDIR/transcripts_all_nonconfunding_region1.txt;

cat $OUTDIR/transcripts_autosomal.txt $OUTDIR/transcripts_hetgam_dropout.txt $OUTDIR/transcripts_heterogametic.txt $OUTDIR/transcripts_partial_hetgam_dropout.txt \
> $OUTDIR/transcripts_all_nonconfunding1.txt;
cat $OUTDIR/all_callable_transcripts.txt | grep -vFf $OUTDIR/transcripts_all_nonconfunding1.txt \
> $OUTDIR/confounding_transcripts1.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable.gtf | grep -vFf $OUTDIR/confounding_transcripts1.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf;

# Remove temporary files
rm $OUTDIR/transcripts_autosomal_temp.txt \
$OUTDIR/transcripts_homogametic_temp.txt \
$OUTDIR/transcripts_hetgam_dropout_temp*.txt \
$OUTDIR/transcripts_heterogametic_temp*.txt;

### Remove transcripts present in more than one strata
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_autosomal.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_autosomal_temp.txt;
cat $BEDS2/$PROJECT\_PAR_3.bed $BEDS2/$PROJECT\_PAR_3_unk.bed \
$BEDS2/$PROJECT\_PAR_5.bed $BEDS2/$PROJECT\_PAR_5_unk.bed \
> $OUTDIR/PAR_regions.bed;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $OUTDIR/PAR_regions.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_PAR_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_PAR_3.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_PAR3_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_PAR_3_unk.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
>> $OUTDIR/transcripts_PAR3unk_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_PAR_5.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
>> $OUTDIR/transcripts_PAR5_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_PAR_5_unk.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
>> $OUTDIR/transcripts_PAR5unk_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_Z.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_Z_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_4A.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_4A_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_3a.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_3a_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_3b.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_3b_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_3c.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_3c_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_5.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_5_temp.txt;
bedtools intersect -a $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf \
-b $BEDS2/$PROJECT\_unknown.bed -f 0.90  | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_unknown_temp.txt;

for i in $(seq 1 1 $(cat $OUTDIR/transcripts_all_nonconfunding_region1.txt | wc -l)); do
  DATA=$(cat $OUTDIR/transcripts_all_nonconfunding_region1.txt | head -n $i | tail -n1);
  FEAT=$(echo $DATA | cut -d ' ' -f1);
  REGION=$(echo $DATA | cut -d ' ' -f2);
  INFO=$(for STRATA in $(echo autosomal PAR Z 4A 3a 3b 3c 5 unknown); do
    echo $(cat $OUTDIR/transcripts_$STRATA\_temp.txt | grep $FEAT) $REGION $STRATA;
  done | awk '{if(NF==3) print}');
  if [ $(echo $INFO | tr ' ' '\n' | wc -l) == 3 ]; then
    if [ $(echo $INFO | cut -d ' ' -f3) == "PAR" ]; then
      PAR=$(for STRATA in $(echo PAR3 PAR3unk PAR5 PAR5unk); do
        cat $OUTDIR/transcripts_$STRATA\_temp.txt | grep $FEAT | awk -v STRATA=$STRATA '{print STRATA}';
      done | awk '{print $NR}');
      INFO=$(echo $(echo $INFO |cut -d' ' -f1,2) $PAR);
    fi;
    echo $INFO | tr ' ' '\t';
  fi;
done | awk -F'\t' '{if($2 == "autosomal" && ($3 == "autosomal" || $3 == "PAR3" || $3 == "PAR3unk" || $3 == "PAR5" || $3 == "PAR5unk")) {print}
else if(($2 == "sex_dropout" || $2 == "sex_phase" || $2 == "partial_dropout") && ($3 == "Z" || $3 == "4A" || $3 == "3a" || $3 == "3b" || $3 == "3c" || $3 == "5")) {print}}' \
> $OUTDIR/transcripts_all_nonconfunding_region2.txt;

cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | cut -f1 \
> $OUTDIR/transcripts_all_nonconfunding2.txt;
cat $OUTDIR/all_callable_transcripts.txt | grep -vFf $OUTDIR/transcripts_all_nonconfunding2.txt \
> $OUTDIR/confounding_transcripts2.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1.gtf | grep -vFf $OUTDIR/confounding_transcripts2.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2.gtf;

# Remove temporary files
rm $OUTDIR/transcripts_autosomal_temp.txt \
$OUTDIR/transcripts_PAR*_temp.txt \
$OUTDIR/transcripts_Z_temp.txt \
$OUTDIR/transcripts_4A_temp.txt \
$OUTDIR/transcripts_3a_temp.txt \
$OUTDIR/transcripts_3b_temp.txt \
$OUTDIR/transcripts_3c_temp.txt \
$OUTDIR/transcripts_5_temp.txt \
$OUTDIR/transcripts_unknown_temp.txt;

## Retain transcripts with partial_dropout over sex_dropout
cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | awk -F'\t' '{if($2 == "sex_dropout") print $1}' | \
awk -F'_' '{print $1"_"$2"_"$3"\t"$1}' \
> $OUTDIR/transcripts_sex_dropout.txt;
cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | awk -F'\t' '{if($2 == "partial_dropout" || $2 == "sex_phase") print $1}' \
> $OUTDIR/transcripts_partial_dropout1.txt;
for i in $(seq 1 1 $(cat $OUTDIR/transcripts_sex_dropout.txt | wc -l)); do
  GENE=$(cat $OUTDIR/transcripts_sex_dropout.txt | head -n$i | tail -n1 | cut -f2);
  PARTIAL=$(cat $OUTDIR/transcripts_partial_dropout1.txt | awk -F'_' -v GENE=$GENE '{if($1==GENE) print}' | wc -l);
  if [ $(echo $PARTIAL) -gt 0 ]; then
    echo $(cat $OUTDIR/transcripts_sex_dropout.txt | head -n$i | tail -n1) remove;
  else
    echo $(cat $OUTDIR/transcripts_sex_dropout.txt | head -n$i | tail -n1) keep;
  fi;
done | awk '{if($3 == "remove") print $1}' \
> $OUTDIR/transcripts_sex_dropout_remove.txt;

## Retain transcripts with sex_phase over partial_dropout
cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | awk -F'\t' '{if($2 == "partial_dropout") print $1}' | \
awk -F'_' '{print $1"_"$2"_"$3"\t"$1}' \
> $OUTDIR/transcripts_partial_dropout2.txt;
cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | awk -F'\t' '{if($2 == "sex_phase") print $1}' \
> $OUTDIR/transcripts_sex_phase.txt;
for i in $(seq 1 1 $(cat $OUTDIR/transcripts_partial_dropout2.txt | wc -l)); do
  GENE=$(cat $OUTDIR/transcripts_partial_dropout2.txt | head -n$i | tail -n1 | cut -f2);
  PARTIAL=$(cat $OUTDIR/transcripts_sex_phase.txt | awk -F'_' -v GENE=$GENE '{if($1==GENE) print}' | wc -l);
  if [ $(echo $PARTIAL) -gt 0 ]; then
    echo $(cat $OUTDIR/transcripts_partial_dropout2.txt | head -n$i | tail -n1) remove;
  else
    echo $(cat $OUTDIR/transcripts_partial_dropout2.txt | head -n$i | tail -n1) keep;
  fi;
done | awk '{if($3 == "remove") print $1}' \
> $OUTDIR/transcripts_partial_dropout_remove.txt;

# Remove transcripts
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2.gtf | \
grep -vFf $OUTDIR/transcripts_sex_dropout_remove.txt | \
grep -vFf $OUTDIR/transcripts_partial_dropout_remove.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf;

## Find transcripts in ingroup species with start codon
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | \
awk -F'\t' '{if($3=="start_codon") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3  | sort -u \
> $OUTDIR/transcripts_start_codon.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt | grep -vFf $OUTDIR/transcripts_start_codon.txt \
> $OUTDIR/transcripts_no_start_codon.txt;

### Transcripts with all CDS lifted over from DB species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | tr -d '"' | tr -d ';' \
> $OUTDIR/features_startstop.txt;
cat $OUTDIR/features_startstop.txt | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_startstop.txt;
cat $OUTDIR/$DB_NAME\_uniqID.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | tr -d '"' | tr -d ';' \
> $OUTDIR/features_$DB_NAME.txt;
cat $OUTDIR/features_$DB_NAME.txt | grep -vFf $OUTDIR/features_startstop.txt | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/transcripts_missing_$DB_NAME.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt | grep -vFf $OUTDIR/transcripts_missing_$DB_NAME.txt \
> $OUTDIR/transcripts_nonmissing.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt  | grep -vFf $OUTDIR/transcripts_nonmissing.txt \
> $OUTDIR/transcripts_missing.txt;

### Transcripts with all CDS lifted over from DB species and with start codon
cat $OUTDIR/transcripts_all_nonconfunding2.txt | grep -vFf $OUTDIR/transcripts_nonmissing.txt | \
grep -vFf $OUTDIR/transcripts_no_start_codon.txt \
> $OUTDIR/transcripts_start_codon_trans_missing.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt | grep -vFf $OUTDIR/transcripts_nonmissing.txt | \
grep -vFf $OUTDIR/transcripts_start_codon_trans_missing.txt \
> $OUTDIR/transcripts_no_start_codon_trans_missing.txt;

# Remove temporary files
rm $OUTDIR/transcripts_start_codon.txt \
$OUTDIR/transcripts_no_start_codon.txt \
$OUTDIR/transcripts_missing_$DB_NAME.txt \
$OUTDIR/transcripts_missing.txt;

## Remove features not overlapping with outgroup species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | awk -F'\t' '{print $NF}' \
> $OUTDIR/loop0.txt;
for i in $(seq 1 1 $(echo $OUTGROUPS | tr ':' '\n' | wc -l)); do
  OUTGROUP=$(echo $OUTGROUPS | tr ':' '\n' | head -n$i | tail -n1 | cut -d' ' -f1);
  cat $(find $OUTDIR -name "loop$(echo $i -1 | bc).txt") | \
  grep -Ff $(find $OUTDIR -name "$OUTGROUP\_features_IDs.tsv") \
  > $OUTDIR/loop$i.txt;
done;
cat $OUTDIR/loop$i.txt \
> $OUTDIR/shared_features.txt;
rm $OUTDIR/loop*.txt;

### Find genes not shared with any outgroup species
cat $OUTDIR/shared_features.txt | awk -F'_' '{if($NF == "CDS") print $1}' | sort -u \
> $OUTDIR/shared_genes.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt | cut -d'_' -f1 | sort -u | \
grep -vxFf $OUTDIR/shared_genes.txt \
> $OUTDIR/nonshared_genes.txt;

### Create annotation file with shared transcripts
cat $OUTDIR/nonshared_genes.txt | awk '{print "gene_id \""$1"\";" }' \
> $OUTDIR/nonshared_genes_temp.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | \
grep -Ff $OUTDIR/nonshared_genes_temp.txt | awk -F'\t' '{print $NF}' | awk -F'_' '{if($NF == "CDS") print}' \
> $OUTDIR/nonshared_features.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | \
grep -vFf $OUTDIR/nonshared_genes_temp.txt | grep -vFf $OUTDIR/shared_features.txt | awk -F'\t' '{print $NF}' \
> $OUTDIR/nonshared_features_remove.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | \
grep -vFf $OUTDIR/nonshared_features_remove.txt | \
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"uniqID \""$10"\";"}' \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf;

# Remove temporary files
rm $OUTDIR/nonshared_genes_temp.txt \
$OUTDIR/nonshared_features_remove.txt;

## Get longest transcripts per gene based on CDS features

### Transcripts with all CDS lifted over from DB species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf | \
grep -vFf $OUTDIR/transcripts_no_start_codon_trans_missing.txt | \
grep -vFf $OUTDIR/transcripts_start_codon_trans_missing.txt | \
awk -F'\t' '{if($3=="CDS") print}' | \
cgat gtf2gtf --method=filter --filter longest-transcript \
-S $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans1.gtf;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans1.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | tr -d '"' | sort -u \
> $OUTDIR/longesttrans1.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans1.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1 | sort -u | awk '{print $1"_"}' \
> $OUTDIR/longesttrans1_genes1.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | \
grep -Ff $OUTDIR/longesttrans1_genes1.txt | tr -d '"' \
> $OUTDIR/longesttrans1_genes2.txt;

### Transcripts with start codon but not all CDS lifted over from DB species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf | \
grep -vFf $OUTDIR/transcripts_no_start_codon_trans_missing.txt | \
grep -vFf $OUTDIR/transcripts_nonmissing.txt | \
grep -vFf $OUTDIR/longesttrans1_genes2.txt | \
awk -F'\t' '{if($3=="CDS") print}' | \
cgat gtf2gtf --method=filter --filter longest-transcript \
-S $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans2.gtf;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans2.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | tr -d '"' | sort -u \
> $OUTDIR/longesttrans2.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans2.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1 | sort -u | awk '{print $1"_"}' \
> $OUTDIR/longesttrans2_genes1.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | \
grep -Ff $OUTDIR/longesttrans2_genes1.txt | tr -d '"' \
> $OUTDIR/longesttrans2_genes2.txt;

### Transcripts without start codon, nor all CDS lifted over from DB species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared.gtf | \
grep -vFf $OUTDIR/transcripts_start_codon_trans_missing.txt | \
grep -vFf $OUTDIR/transcripts_nonmissing.txt | \
grep -vFf $OUTDIR/longesttrans1_genes2.txt | \
grep -vFf $OUTDIR/longesttrans2_genes2.txt | \
awk -F'\t' '{if($3=="CDS") print}' | \
cgat gtf2gtf --method=filter --filter longest-transcript \
-S $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans3.gtf;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans3.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | tr -d '"' | sort -u \
> $OUTDIR/longesttrans3.txt;

### Merge annotation files
cat $OUTDIR/longesttrans1.txt $OUTDIR/longesttrans2.txt $OUTDIR/longesttrans3.txt | sort \
> $OUTDIR/shared_longesttrans.txt;
cat $OUTDIR/transcripts_all_nonconfunding2.txt | grep -vFf $OUTDIR/shared_longesttrans.txt \
> $OUTDIR/nonlongesttrans.txt;
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region.gtf | \
grep -vFf $OUTDIR/nonlongesttrans.txt \
> $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans.gtf;

# Remove temporary files
rm $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans1.gtf \
$OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans2.gtf \
$OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans3.gtf \
$OUTDIR/longesttrans1_genes1.txt \
$OUTDIR/longesttrans1_genes2.txt \
$OUTDIR/longesttrans2_genes1.txt \
$OUTDIR/longesttrans2_genes2.txt;

## Find transcripts with corresponding features in database species
cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans.gtf | \
awk -F'\t' '{if($3=="CDS") print}' | awk '{print $NF}' | cut -d'_' -f1,2,3 | sort -u \
> $OUTDIR/longest_transcripts.txt;
cat $OUTDIR/$DB_NAME\_uniqID.gtf | awk -F'\t' '{if($3=="CDS") print}' | \
cut -f1,4,5,7,10 | grep -Ff $OUTDIR/longest_transcripts.txt \
> $OUTDIR/$DB_NAME\_features.tsv;

## Collect corresponding data in ingroup species
for i in $(seq 1 1 $(cat $OUTDIR/$DB_NAME\_features.tsv | wc -l)); do
  DB=$(cat $OUTDIR/$DB_NAME\_features.tsv | head -n$i | tail -n1);
  FEAT=$(echo $DB | cut -d' ' -f 5);
  DATA=$(cat $OUTDIR/$SPEC_NAME\_uniqID_uniqfeats_startstop_callable_nonconf1_nonconf2_region_shared_longesttrans.gtf | \
  grep $FEAT | cut -f 1,4,5,7);
  TRANS=$(echo $FEAT | cut -d'_' -f2,3);
  REGION=$(cat $OUTDIR/transcripts_all_nonconfunding_region2.txt | grep $TRANS | cut -f2,3);
  if [ $(echo $REGION | cut -d' ' -f1) == "partial_dropout" ]; then
    CALLABLE=$(cat $OUTDIR/transcripts_partial_dropout_feature_dropout.txt | grep $FEAT);
    if [ $(echo $CALLABLE | tr -d '\n' | wc -c) -gt 0 ]; then
      CALLABLE="partial_dropout";
    else
      CALLABLE="callable";
    fi;
  else
    CALLABLE="callable";
  fi;
  SHARED=$(cat $OUTDIR/nonshared_features.txt | grep $FEAT);
  if [ $(echo $SHARED | tr -d '\n' | wc -c) -gt 0 ]; then
    SHARED="nonshared";
  else
    SHARED="shared";
  fi;
  if [ $(echo $DATA | tr -d '\n' | wc -c) -gt 0 ]; then
    echo -e "${DB}\t${DATA}\t${CALLABLE}\t${SHARED}\t${REGION}\t${FEAT}";
  else
    echo -e "${DB}\tmissing\tmissing\tmissing\tmissing\tmissing\tmissing\t${REGION}\t${FEAT}";
  fi;
done | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' | \
grep -v "unassigned_transcript" \
> $OUTDIR/$PROJECT\_features.tsv;
