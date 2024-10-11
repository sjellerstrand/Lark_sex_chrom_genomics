#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J ancestralZ_strata
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
GENES=$WORKDIR/D4_align_shared_genes/metadata/genes_info.tsv;
CHICKEN_PROT=$MAINDIR/data/reference/Gallus_gallus/GCF_000002315.4_Gallus_gallus-5.0_protein.faa;
CHICKEN_GTF=$MAINDIR/data/reference/Gallus_gallus/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff;
FLYCATCHER_PROT=$MAINDIR/data/reference/Ficedula_albicollis/GCF_000247815.1_FicAlb1.5_protein.faa;
FLYCATCHER_GTF=$MAINDIR/data/reference/Ficedula_albicollis/GCF_000247815.1_FicAlb1.5_genomic.gtf;
DB_NAME=Tgut;
DB_REF=$MAINDIR/data/Skylark_2021/B3_annotation_lift_over/$DB_NAME.fa;
DB_TRANS=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_rna.fna;
DB_GTF=$MAINDIR/data/reference/Taeniopygia_guttata/GCF_003957565.2_bTaeGut1.4.pri_genomic.no_W.gtf;
OSTRICH_GTF=$WORKDIR/B3_annotation_lift_over_Scam/Scam_Tgut_liftover.gtf_polished;
EMU_GTF=$WORKDIR/B3_annotation_lift_over_Dnov/Dnov_Tgut_liftover.gtf_polished;
METADATA=$WORKDIR/metadata;
ANCESTRAL_STRATA=$METADATA/Ancestral_strata_genes.tsv;
OSTRICH_Z_STRUCTURE=$METADATA/Ostrich_Z_structure.tsv;

### Load modules
module load bioinfo-tools blast/2.15.0+;

### Create folders
mkdir $OUTDIR/ancestralZ_strata;
OUTDIR=$OUTDIR/ancestralZ_strata;
mkdir $OUTDIR/database \
$OUTDIR/matches \
$OUTDIR/metadata;

cat $GENES | awk -F'\t' '{if($5 == "Z" || $5 == "Z") print $1"\t"$2}' \
> $OUTDIR/metadata/Ancestral_Z_genes.tsv;

echo -e "Gene\tTranscript\tXu_et_al_2019_gene\tTgut_stratum\tBLAST_Ggul_gene\tBlast_Ggul_ID\tBLAST_Ggul_stratum\tBLAST_Ggul_match_attempt\tBLAST_Falc_gene\tBLAST_Falc_stratum\tBLAST_Falc_match_attempt\tScam_scaffold\tScam_start\tScam_end\tDnov_scaffold\tDnov_start\tDnov_end" \
> $OUTDIR/metadata/strata_search.tsv;

# Create blast database from Chicken
cp $CHICKEN_PROT $OUTDIR/database/chicken_protein.faa;
cp $DB_TRANS $OUTDIR/database/$DB_NAME\_rna.fna;
makeblastdb -in $OUTDIR/database/chicken_protein.faa -dbtype prot;

# Create blast database from Collared flycatcher
cp $FLYCATCHER_PROT $OUTDIR/database/flycatcher_protein.faa;
cp $DB_TRANS $OUTDIR/database/$DB_NAME\_rna.fna;
makeblastdb -in $OUTDIR/database/flycatcher_protein.faa -dbtype prot;

# Match genes to ancestral strata
for i in $(seq 1 1 $(cat $OUTDIR/metadata/Ancestral_Z_genes.tsv | wc -l)); do
  GENE=$(cat $OUTDIR/metadata/Ancestral_Z_genes.tsv | head -n$i | tail -n1 | cut -f1);
  TRANS=$(cat $OUTDIR/metadata/Ancestral_Z_genes.tsv | head -n$i | tail -n1 | cut -f2);
  mkdir $OUTDIR/matches/$GENE\_$TRANS;

  # Match DB name to Xu et al. 2019 gene list
  ANC_MATCH=$(cat $ANCESTRAL_STRATA | awk -F'\t' -v GENE=$GENE '{if(tolower($1) == tolower(GENE)) print $1"\t"$2}');
  if [ $(echo $ANC_MATCH | tr -d '\n' | wc -c) -eq 0 ]; then
    ANC_MATCH=$(echo -e "NA\tNA");
  fi;

  # Extract database sequence
  cat $OUTDIR/database/$DB_NAME\_rna.fna | awk -v TRANS=$TRANS 'BEGIN {printline="F"} \
  {if($0 ~ "^" ">"TRANS) {printline="T"} else if($0 ~ "^" ">") {printline="F"}; if(printline=="T") {print}}' \
  > $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS.fasta;

  # Blastx against Chicken
  blastx -query $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS.fasta -db $OUTDIR/database/chicken_protein.faa \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle" \
  -out $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_chicken_match -evalue 1e-6 -max_target_seqs 10;

  # Evaluate top match. Try at least 10 attempts
  if [ $(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_chicken_match | wc -l) -gt 10 ]; then
    ATTEMPTS=10;
  else
    ATTEMPTS=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_chicken_match | wc -l);
  fi;
  if [ $(echo $ATTEMPTS) -eq 0 ]; then
    BLAST_GGUL_MATCH=$(echo -e "NA\tNA\tNA\t0");
  else
    for ATTEMPT in $(seq 1 1 $ATTEMPTS); do
      MATCH=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_chicken_match | head -n$ATTEMPT | tail -n1 | cut -f2);
      GENE_MATCH=$(cat $CHICKEN_GTF | grep $MATCH | awk -F'\t' '{if($3 == "CDS") print $9}' | \
      awk -F'[=;]' -v MATCH=$MATCH '{for(i=1; i<=NF; i++) {if($i=="protein_id") {FIELD=i+1; break}}; if($FIELD==MATCH) {print; exit}}' | \
      awk -F'[=;]' '{for(i=1; i<=NF; i++) {if($i=="Dbxref") {ID=i+1}; if($i=="gene") {GENE=i+1}}; print $GENE"\t"$ID}' | awk -F'[\t:,]' '{print $1"\t"$3}');
      if [ $(echo $GENE_MATCH | tr -d '\n' | wc -c) -gt 0 ]; then
        GENE_ID=$(echo $GENE_MATCH | cut -d' ' -f2);
        STRATUM=$(cat $ANCESTRAL_STRATA | awk -F'\t' -v GENE_ID=$GENE_ID '{if(tolower($3) == tolower(GENE_ID)) print $2}');
        if [ $(echo $STRATUM | tr -d '\n' | wc -c) -eq 0 ]; then
          STRATUM=$(echo "NA");
        fi;
        BLAST_GGUL_MATCH=$(echo -e "$GENE_MATCH\t$STRATUM\t$ATTEMPT");
        break;
      fi;
    done;
    if [ $(echo $GENE_MATCH | tr -d '\n' | wc -c) -eq 0 ]; then
      BLAST_GGUL_MATCH=$(echo -e "NA\tNA\tNA\t$ATTEMPT");
    fi;
  fi;

  # Blastx against Collared flycatcher
  blastx -query $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS.fasta -db $OUTDIR/database/flycatcher_protein.faa \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen evalue bitscore stitle" \
  -out $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_flycatcher_match -evalue 1e-6 -max_target_seqs 10;

  # Evaluate top match. Try at least 10 attempts
  if [ $(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_flycatcher_match | wc -l) -gt 10 ]; then
    ATTEMPTS=10;
  else
    ATTEMPTS=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_flycatcher_match | wc -l);
  fi;
  if [ $(echo $ATTEMPTS) -eq 0 ]; then
    BLAST_FALC_MATCH=$(echo -e "NA\tNA\t0");
  else
    for ATTEMPT in $(seq 1 1 $ATTEMPTS); do
      MATCH=$(cat $OUTDIR/matches/$GENE\_$TRANS/$GENE\_$TRANS\_flycatcher_match | head -n$ATTEMPT | tail -n1 | cut -f2);
      GENE_MATCH=$(cat $FLYCATCHER_GTF | grep $MATCH | awk -F'\t' '{if($3 == "CDS") print $9}' | \
      awk -F'"' -v MATCH=$MATCH '{for(i=1; i<=NF; i++) {if($i=="; protein_id ") {FIELD=i+1; break}}; if($FIELD==MATCH) {print $2; exit}}' | head -n1);
      if [ $(echo $GENE_MATCH | tr -d '\n' | wc -c) -gt 0 ]; then
        STRATUM=$(cat $ANCESTRAL_STRATA | awk -F'\t' -v GENE_MATCH=$GENE_MATCH '{if(tolower($1) == tolower(GENE_MATCH)) print $2}');
        if [ $(echo $STRATUM | tr -d '\n' | wc -c) -eq 0 ]; then
          STRATUM=$(echo "NA");
        fi;
        BLAST_FALC_MATCH=$(echo -e "$GENE_MATCH\t$STRATUM\t$ATTEMPT");
        break;
      fi;
    done;
    if [ $(echo $GENE_MATCH | tr -d '\n' | wc -c) -eq 0 ]; then
      BLAST_FALC_MATCH=$(echo -e "NA\tNA\t$ATTEMPT");
    fi;
  fi;

  # Find corresponding position in Ostrich
  OSTRICH_POS=$(cat $OSTRICH_GTF | grep $TRANS | awk -F'[\t"]' -v TRANS=$TRANS '{if($3 == "transcript" && $12 == TRANS) print $1"\t"$4"\t"$5}');
  if [ $(echo $OSTRICH_POS | tr -d '\n' | wc -c) -eq 0 ]; then
    OSTRICH_POS=$(echo -e "NA\tNA\tNA");
  else
    Z_POS=$(cat $OSTRICH_Z_STRUCTURE | grep $(echo $OSTRICH_POS | cut -d' ' -f1));
    if [ $(echo $Z_POS | tr -d '\n' | wc -c) -eq 0 ]; then
      OSTRICH_POS=$(echo -e "NA\tNA\tNA");
    else
      N_ROWS=$(echo $Z_POS | cut -d' ' -f1);
      SUM_POS=$(cat $OSTRICH_Z_STRUCTURE | tail -n+2 | head -n$N_ROWS | awk -F'\t' -v N_ROWS=$N_ROWS 'BEGIN {sum=0} {if(NR == N_ROWS && $5=="+") {sum=sum} else {sum+=$4}} END {print sum}');
      if [ $(echo $Z_POS | cut -d' ' -f5) == "+" ]; then
        OSTRICH_POS=$(echo $OSTRICH_POS | awk -v SUM_POS=$SUM_POS '{print "Z\t"SUM_POS+$2"\t"SUM_POS+$3}');
      else
        OSTRICH_POS=$(echo $OSTRICH_POS | awk -v SUM_POS=$SUM_POS '{print "Z\t"SUM_POS-$2"\t"SUM_POS-$3}');
      fi;
    fi;
  fi;

  # Find corresponding position in Emu
  EMU_POS=$(cat $EMU_GTF | grep $TRANS | awk -F'[\t"]' -v TRANS=$TRANS '{if($3 == "transcript" && $12 == TRANS) print $1"\t"$4"\t"$5}');
  if [ $(echo $EMU_POS | tr -d '\n' | wc -c) -eq 0 ]; then
    EMU_POS=$(echo -e "NA\tNA\tNA");
  fi;

  echo -e "$GENE\t$TRANS\t$ANC_MATCH\t$BLAST_GGUL_MATCH\t$BLAST_FALC_MATCH\t$OSTRICH_POS\t$EMU_POS" | tr -s ' ' | tr ' ' '\t';

done >> $OUTDIR/metadata/strata_search.tsv;

# Decide ancestral strata
echo -e "Gene\tTranscript\tMatched_gene\tStratum\tConfidence_rank" \
> $OUTDIR/metadata/strata_match.tsv;
cat $OUTDIR/metadata/strata_search.tsv | tail -n+2 | \
awk -F'\t' '{MATCH="NA"; STRATUM="NA"; RANK="NA";
if($7!="NA" && $8 == 1) {MATCH=$5; STRATUM=$7; RANK=1} \
else if($4!="NA") {MATCH=$3; STRATUM=$4; RANK=2} \
else if($10!="NA" && $11 == 1) {MATCH=$9; STRATUM=$10; RANK=3} \
else if($12=="Z" && $15 == "Z") {MATCH="NA"; RANK=4;
if($16 <= 24118484 && $17 <= 24118484) {STRATUM="S0"} \
else if($16 >= 24238348 && $17 >= 24238348 && $16 <= 55994283 && $17 <= 55994283) {STRATUM="S1"} \
else if($16 >= 56117152 && $17 >= 56117152 && $16 <= 75469006 && $17 <= 75469006) {STRATUM="S2"} \
else if($16 >= 75514489 && $17 >= 75514489) {STRATUM="S3"}} \
print $1"\t"$2"\t"MATCH"\t"STRATUM"\t"RANK}' \
>> $OUTDIR/metadata/strata_match.tsv;
