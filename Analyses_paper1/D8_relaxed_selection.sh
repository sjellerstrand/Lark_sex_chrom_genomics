#!/bin/bash -l

#SBATCH -A naiss2024-5-340
#SBATCH -p shared
#SBATCH -n 2
#SBATCH -t 20:00:00
#SBATCH -J relaxed_selection
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/cfs/klemming/projects/supr/snic2020-2-25/user_data/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/$PROJECT1;
GENES=$WORKDIR1/organise_data1/$PROJECT1\_$PROJECT2\_organised_data1.tsv;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
TREE=$METADATA1/tree_topologies/Skylark_2021_Rasolark_2021_sex_phase_sex_phase_females.nwk;

DB_NAME=Tgut;
DB_NAME=Tgut;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");

### Activate conda environment
conda activate hyphy;

### Create folders
mkdir $OUTDIR/D8_relaxed_selection;
OUTDIR=$OUTDIR/D8_relaxed_selection;
mkdir $OUTDIR/metadata \
$OUTDIR/genes;

# Subset genes
cat $GENES | awk -F'\t' '{if($98 == "OK" && $100=="OK" && $101=="OK" && $102=="OK" && \
$4!="S0" && $4!="S1" && $4!="S2" && $4!="S3" && $4!="Ancestral unknown" && $3=="sex_phase") print $1}' | \
sort | uniq -D | uniq \
> $OUTDIR/metadata/genes_in.tsv;
GENES2=$OUTDIR/metadata/genes_in.tsv;

# Annotate branches in tree
cat $TREE | sed -E 's/\b(BranchW)\b/&{Test}/g' | \
sed -E 's/\b(BranchZ)\b/&{Reference}/g' \
> $OUTDIR/metadata/tree.nwk;
TREE=$OUTDIR/metadata/tree.nwk;

for GENE in $(cat $GENES2); do

  echo $GENE;
  mkdir $OUTDIR/genes/$GENE;
  cd $OUTDIR/genes/$GENE;
  
  hyphy relax \
  --alignment $WORKDIR1/D4_align_shared_genes/shared_genes/sex_phase_sex_phase/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_ZW_ZW.fasta \
  --code Universal --tree $TREE --models Minimal --reference Reference --test Test \
  --output $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.json \
  --save-fit $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.txt;


done;

# Retrieve output data
echo -e "geneID\tintensification_parameter\tp-value\tcorrected_p-value" \
> $OUTDIR/metadata/relax_results.tsv;

for GENE in $(cat $GENES2); do

  INTPAR=$(echo $(jq '.["test results"]["relaxation or intensification parameter"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.json));
  echo -e "${GENE}\t${INTPAR}";


echo $(jq '.["branch attributes"]["0"]["BranchZ"]["Uncorrected P-value"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.json)

done >> $OUTDIR/metadata/relax_results.tsv;

### Correct for multiple tests. (Not all tests are actually performed).
NGENES=$(cat $OUTDIR/metadata/relax_results.tsv | tail -n+2| awk -F'\t' 'BEGIN {Ngen=0} {if($2 != "") Ngen+=1} END {print Ngen}');
echo -e "geneID\tintensification_parameter\tp-value\tcorrected_p-value" \
> $OUTDIR/metadata/relax_results.tsv;

for GENE in $(cat $GENES2); do

  INTPAR=$(echo $(jq '.["test results"]["relaxation or intensification parameter"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.json));
  PVAL=$(echo $(jq '.["test results"]["p-value"]' \
  $OUTDIR/genes/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_RELAX.json));
  PVALCORR=$(echo $PVAL | awk -v NGENES=$NGENES '{pcorr=$1*NGENES; if(pcorr > 1) {pcorr=1}; print pcorr}');
  echo -e "${GENE}\t${INTPAR}\t${PVAL}\t${PVALCORR}";

done >> $OUTDIR/metadata/relax_results.tsv;
