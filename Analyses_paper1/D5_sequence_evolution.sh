#!/bin/bash -l

#SBATCH -A naiss2024-5-92
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J sequence_evolution
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT1=Skylark_2021;
PROJECT2=Rasolark_2021;
WORKDIR1=$MAINDIR/data/$PROJECT1;
WORKDIR2=$MAINDIR/data/$PROJECT2;
OUTDIR=$MAINDIR/working/Simon/$PROJECT1;
GENES=$WORKDIR1/D4_align_shared_genes/metadata/genes_info.tsv;
DB_NAME=Tgut;
METADATA1=$WORKDIR1/metadata;
METADATA2=$WORKDIR2/metadata;
DB_NAME=Tgut;
OUTGROUPS=$(echo -e \
"Pmaj \
${MAINDIR}/data/reference/Parus_major/GCF_001522545.3_Parus_major1.1_genomic.fasta \
${WORKDIR}/D1_filter_genes/Pmaj_uniqID.gtf
:\
Falb \
${MAINDIR}/data/reference/Ficedula_albicollis/GCA_000247815.2_FicAlb1.5_genomic.fasta \
${WORKDIR}/D1_filter_genes/Falb_uniqID.gtf");
hyphy_analyses=/crex/proj/snic2020-2-25/bin/hyphy-analyses;

### Load modules
module load bioinfo-tools SeqKit/0.15.0 gnuparallel/20180822;

### Activate conda environment
conda activate hyphy;

### Create folders
mkdir $OUTDIR/D5_sequence_evolution;
OUTDIR=$OUTDIR/D5_sequence_evolution;
mkdir $OUTDIR/metadata \
$OUTDIR/shared_genes \
$OUTDIR/shared_genes/autosomal_autosomal \
$OUTDIR/shared_genes/sex_dropout_sex_dropout \
$OUTDIR/shared_genes/sex_phase_sex_phase \
$OUTDIR/shared_genes/sex_phase_sex_dropout \
$OUTDIR/shared_genes/sex_dropout_sex_phase \
$OUTDIR/shared_genes/sex_phase_autosomal;

# Set up metadata files
OUTGROUPS=$(echo $OUTGROUPS | tr ':' '\n' | cut -d' ' -f1 | tr '\n' '|');
OUTGROUPS=${OUTGROUPS::-1};
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_A_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT1\_female_A_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT1\_female_Z_sampleID.txt;
cat $WORKDIR1/D3_codon_aware_alignment/metadata/female_W_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT1\_female_W_sampleID.txt;
cat $WORKDIR2/D3_codon_aware_alignment/metadata/female_A_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT2\_female_A_sampleID.txt;
cat $WORKDIR2/D3_codon_aware_alignment/metadata/female_Z_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT2\_female_Z_sampleID.txt;
cat $WORKDIR2/D3_codon_aware_alignment/metadata/female_W_sampleID.txt | \
awk -v DB_NAME=$DB_NAME -v OUTGROUPS=$OUTGROUPS '{if($0 !~ DB_NAME && $0 !~ OUTGROUPS) print}' \
>> $OUTDIR/metadata/$PROJECT2\_female_W_sampleID.txt;

export WORKDIR1 WORKDIR2 OUTDIR FEATURES1 FEATURES2 PROJECT1 PROJECT2 METADATA1 METADATA2 hyphy_analyses;

sequence_evolution() {

  GENE=$1;
  TRANS=$2;
  REGION0_1=$3;
  REGION0_2=$4;
  STRATA1=$5;
  STRATA2=$6;

  # If gene is not shared between species, abort
  if [ ${14} != $PROJECT1\_$PROJECT2 ]; then
    echo $GENE not shared between species;
    exit;
  elif [ ${10} == 0 ]; then
    echo $GENE trimmed length is 0;
    exit;
  fi;

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

  echo $GENE;
  mkdir $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE;

  # Choose branches of interest
  if [ $REGION1 == "autosomal" ] && [ $REGION2 == "autosomal" ]; then
      BRANCHES=$(echo BranchA BranchSkyA BranchRasoA);
      SUFFIX="A_A";
      GROUPS1=$(echo A);
      GROUPS2=$(echo A);
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="Z_Z";
    GROUPS1=$(echo Z);
    GROUPS2=$(echo Z);
    if [ $STRATA1 == "Z" ] && [ $STRATA2 == "Z" ]; then
      BRANCHES=$(echo BranchAZ BranchSkyZ BranchRasoZ);
    else
      BRANCHES=$(echo BranchAZ BranchSkyZ BranchRasoZ);
    fi;
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="ZW_ZW";
    GROUPS1=$(echo Z W);
    GROUPS2=$(echo Z W);
    if [ $STRATA1 == "Z" ] && [ $STRATA2 == "Z" ]; then
      BRANCHES=$(echo BranchZ BranchSkyZ BranchRasoZ BranchW BranchSkyW BranchRasoW);
    else
      BRANCHES=$(echo BranchA BranchZ BranchSkyZ BranchRasoZ BranchW BranchSkyW BranchRasoW);
    fi;
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "sex_dropout" ]; then
    SUFFIX="ZW_Z";
    GROUPS1=$(echo Z W);
    GROUPS2=$(echo Z);
    if [ $STRATA1 == "Z" ] && [ $STRATA2 == "Z" ]; then
      BRANCHES=$(echo BranchZ BranchSkyZ BranchRasoZ BranchSkyW);
    else
      BRANCHES=$(echo BranchA BranchZ BranchSkyZ BranchRasoZ BranchSkyW)
    fi;
  elif [ $REGION1 == "sex_dropout" ] && [ $REGION2 == "sex_phase" ]; then
    SUFFIX="Z_ZW";
    GROUPS1=$(echo Z);
    GROUPS2=$(echo Z W);
    if [ $STRATA1 == "Z" ] && [ $STRATA2 == "Z" ]; then
      BRANCHES=$(echo BranchZ BranchSkyZ BranchRasoZ BranchRasoW);
    else
      BRANCHES=$(echo BranchA BranchZ BranchSkyZ BranchRasoZ BranchRasoW)
    fi;
  elif [ $REGION1 == "sex_phase" ] && [ $REGION2 == "autosomal" ]; then
    SUFFIX="ZW_A";
    GROUPS1=$(echo Z W);
    GROUPS2=$(echo A);
    BRANCHES=$(echo BranchA BranchRasoA BranchSkyZ BranchSkyW BranchO Node1);
  fi;

  # Choose tree topology
  if [ $STRATA1 == "Z" ] && [ $STRATA2 == "Z" ]; then
    TREE=$METADATA1/tree_topologies/$PROJECT1\_$PROJECT2\_$REGION1\_$REGION2\_ancestralZ_females.nwk;
  else
    TREE=$METADATA1/tree_topologies/$PROJECT1\_$PROJECT2\_$REGION1\_$REGION2\_females.nwk;
  fi;

  # Fit codon substitution model to alignment
  cat $WORKDIR1/D4_align_shared_genes/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females_$SUFFIX.fasta \
  > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females.fasta;
  hyphy $hyphy_analyses/FitMG94/FitMG94.bf "ENV=TOLERATE_NUMERICAL_ERRORS=1;" \
  --alignment $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females.fasta \
  --rooted Yes --code Universal --tree $TREE --type local --lrt No --frequencies CF3x4 \
  --output $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json \
  --save-fit $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_fit.txt;

  # Retrieve within group pNpS
  touch $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt;

  # Project1
  for GROUP in $(echo $GROUPS1); do

    # Subset alignment to group sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females.fasta | \
    seqkit grep -f $OUTDIR/metadata/$PROJECT1\_female_$GROUP\_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$GENE\_alignment_NT_females_$GROUP\_ingroup.fasta;

    # Fit codon substitution model to alignment
    hyphy $hyphy_analyses/FitMG94/FitMG94.bf "ENV=TOLERATE_NUMERICAL_ERRORS=1;" \
    --alignment $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$GENE\_alignment_NT_females_$GROUP\_ingroup.fasta \
    --code Universal --tree $METADATA1/tree_topologies/$PROJECT1\_$REGION1\_females_ingroup_$GROUP.nwk --type global --lrt No --frequencies CF3x4 \
    --output $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$GENE\_MG94_$GROUP\_ingroup.json \
    --save-fit $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$GENE\_MG94_fit_$GROUP\_ingroup.txt;

    # Extract omega
    OMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$GENE\_MG94_$GROUP\_ingroup.json | \
    grep "non-synonymous/synonymous rate ratio" | tail -n1 | cut -d':' -f2 | cut -d',' -f1);

    # Report data
    echo -e "${PROJECT1}\t${GROUP}\t${OMEGA}" \
    >> $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt;

  done;

  # Project2
  for GROUP in $(echo $GROUPS2); do

    # Subset alignment to group sequences only
    cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_alignment_NT_females.fasta | \
    seqkit grep -f $OUTDIR/metadata/$PROJECT2\_female_$GROUP\_sampleID.txt | tr -d ' ' | \
    awk '{if($0 ~ /^>/) printf("\n"$0"\n"); else printf($0)} END {printf("\n")}' | tail -n+2 \
    > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT2\_$GENE\_alignment_NT_females_$GROUP\_ingroup.fasta;

    # Fit codon substitution model to alignment
    hyphy $hyphy_analyses/FitMG94/FitMG94.bf "ENV=TOLERATE_NUMERICAL_ERRORS=1;" \
    --alignment $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT2\_$GENE\_alignment_NT_females_$GROUP\_ingroup.fasta \
    --code Universal --tree $METADATA2/tree_topologies/$PROJECT2\_$REGION2\_females_ingroup_$GROUP.nwk --type global --lrt No --frequencies CF3x4 \
    --output $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT2\_$GENE\_MG94_$GROUP\_ingroup.json \
    --save-fit $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT2\_$GENE\_MG94_fit_$GROUP\_ingroup.txt;

    # Extract omega
    OMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT2\_$GENE\_MG94_$GROUP\_ingroup.json | \
    grep "non-synonymous/synonymous rate ratio" | tail -n1 | cut -d':' -f2 | cut -d',' -f1);

    # Report data
    echo -e "${PROJECT2}\t${GROUP}\t${OMEGA}" \
    >> $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt;

  done;

  # Extract branch specific data
  for BRANCH in $BRANCHES; do

    if [ ${#BRANCH} -gt 8 ]; then
      if [ $(echo $BRANCH | grep Sky | wc -l) -gt 0 ]; then
        if [ $(echo $BRANCH | rev | grep ^A | wc -l) -gt 0 ]; then
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt | awk -F'\t' '{if($1 ~ /^Sky/ && $2=="A") print $3}');
        elif [ $(echo $BRANCH | rev | grep ^Z | wc -l) -gt 0 ]; then
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt | awk -F'\t' '{if($1 ~ /^Sky/ && $2=="Z") print $3}');
        else
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt | awk -F'\t' '{if($1 ~ /^Sky/ && $2=="W") print $3}');
        fi;
      elif [ $(echo $BRANCH | grep Raso | wc -l) -gt 0 ]; then
        if [ $(echo $BRANCH | rev | grep ^A | wc -l) -gt 0 ]; then
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt | awk -F'\t' '{if($1 ~ /^Raso/ && $2=="A") print $3}');
        elif [ $(echo $BRANCH | rev | grep ^Z | wc -l) -gt 0 ]; then
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt |  awk -F'\t' '{if($1 ~ /^Raso/ && $2=="Z") print $3}');
        else
          INGROUPOMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt |  awk -F'\t' '{if($1 ~ /^Raso/ && $2=="W") print $3}');
        fi;
      fi;
      if [ -z "${INGROUPOMEGA}" ]; then
        INGROUPOMEGA="NA";
      fi;
    else
      INGROUPOMEGA="NA";
    fi;

    if [ $(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | grep $BRANCH | wc -l) == 0 ]; then

      # Report zero-length branch data as zeroes and NA
      echo -e "${GENE}\t${BRANCH}\t0\t0\t0\tNA\t0\t0\tNA\t$INGROUPOMEGA";

    else

      # Extract non-synonymous substitutions per codon
      NONSYN=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | \
      tr -d '[\n ]' | tr '[{}]' '\n' | grep -E ,\"$BRANCH\": -A3 | \
      head -n4 | tail -n1 | cut -d':' -f6 | cut -d',' -f1);

      # Extract synonymous substitutions per codon
      SYN=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | \
      tr -d '[\n ]' | tr '[{}]' '\n' | grep -E ,\"$BRANCH\": -A3 | \
      head -n4 | tail -n1 | cut -d':' -f7 | cut -d',' -f1);

      # Extract branch length
      LENGTH=$(echo $NONSYN $SYN | awk '{print $1+$2}');

      # Extract dN
      dN=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | \
      tr -d '[\n ]' | tr '[{}]' '\n' | grep -E ,\"$BRANCH\": -A3 | \
      head -n4 | tail -n1 | cut -d':' -f4 | cut -d',' -f1);

      # Extract dS
      dS=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | \
      tr -d '[\n ]' | tr '[{}]' '\n' | grep -E ,\"$BRANCH\": -A3 | \
      head -n4 | tail -n1 | cut -d':' -f5 | cut -d',' -f1);

      # Calculate dN/dS
      if [ $dS == 0 ]; then
        dNdS=$(echo NA);
      else
        dNdS=$(echo $dN $dS | awk '{print $1/$2}');
      fi;

      # Extract beta
      BETA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_fit.txt | \
      grep $BRANCH.beta | cut -d'=' -f2 | cut -d ';' -f1);

      # Extract alpha
      ALPHA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_fit.txt | \
      grep $BRANCH.alpha | cut -d'=' -f2 | cut -d ';' -f1);

      # Extract omega
      OMEGA=$(cat $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94.json | \
      tr -d '[\n ]' | tr '[{}]' '\n' | grep -E ,\"$BRANCH\": -A3 | \
      head -n3 | tail -n1 | cut -d':' -f3 | cut -d',' -f1);

      # Report data
      echo -e "${GENE}\t${BRANCH}\t${LENGTH}\t${dN}\t${dS}\t${dNdS}\t${BETA}\t${ALPHA}\t${OMEGA}\t${INGROUPOMEGA}";

    fi;

  done > $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data.txt;

  # Remove temporary files
  rm -r  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/*fasta \
  $OUTDIR/shared_genes/$REGION1\_$REGION2/$GENE/$PROJECT1\_$PROJECT2\_$GENE\_MG94_data_ingroup.txt;

};

## Excecute function in parallell
export -f sequence_evolution;
parallel --colsep '\t' 'sequence_evolution {}' :::: $GENES;

# Make file with info of each gene alignment
cat $(find $OUTDIR/shared_genes -name "*_MG94_data.txt") | sort -k1 \
> $OUTDIR/metadata/genes_sequence_evolution_info.tsv;
rm $(find $OUTDIR/shared_genes -name "*_MG94_data.txt");

