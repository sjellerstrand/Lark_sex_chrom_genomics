#!/bin/bash -l

#SBATCH -A snic2021-5-469
#SBATCH -p node
#SBATCH -n 1
#SBATCH -C mem256GB
#SBATCH -t 06:00:00
#SBATCH -J prepare_reference
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-2-25/nobackup/simon/Sylvioidea;
PROJECT=Skylark_2021;
WORKDIR=$MAINDIR/data/$PROJECT;
OUTDIR=$MAINDIR/working/Simon/$PROJECT;
REF=$MAINDIR/data/reference/Alauda_arvensis\
/GCA_902810485.1_skylark_genome_genomic.fasta;
REPEAT_QUERY=Aves;
convert_IUPAC=/crex/proj/snic2020-2-25/bin\
/convert_ambiguity-accessed-2021-12-20/convert_ambiguity_SJE.py;

### Load modules
module load bioinfo-tools biopython/1.76-py3 \
samtools/1.14 bwa/0.7.17 RepeatMasker/4.1.0;

### Create folders
mkdir $OUTDIR/B1_prepare_reference;
OUTDIR=$OUTDIR/B1_prepare_reference;
mkdir $OUTDIR/repeatmasker;

### Check for IUPAC ambiguity codes.
IUPAC=$(cat $REF | grep -v "^>" | \
grep [KkMmRrYySsWwVvHhDdBb] | awk 'BEGIN {IUPAC=0} {if($1 > 1) IUPAC=1; exit} END {print IUPAC}');

### If IUPAC ambiguity codes, Pick random base for K, M, R, Y,S and W. Set V, H, D and B to N.
if [ $IUPAC -gt 0 ]; then
    echo "IUPAC ambiguity codes found. Reference will be modified"
    python $convert_IUPAC $REF $OUTDIR/$(echo $REF | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)\_no_IUPAC.fasta;
    REF=$OUTDIR/$(echo $REF | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)\_no_IUPAC.fasta;
    samtools faidx $REF;
    bwa index $REF;
else
    echo "No IUPAC ambiguity codes. Reference will not be modified"
    cp $REF $OUTDIR/$(echo $REF | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)\_no_IUPAC.fasta;
    REF=$OUTDIR/$(echo $REF | rev | cut -d'/' -f1 | cut -d'.' -f2- | rev)\_no_IUPAC.fasta;
    samtools faidx $REF;
    bwa index $REF;
fi;

### Mask repeats
RepeatMasker -par 20 -species $REPEAT_QUERY -dir $OUTDIR/repeatmasker $REF;
cat $OUTDIR/repeatmasker/$(echo $REF | rev | cut -d'/' -f1 | rev).out | \
tr -s ' ' | tail -n+4 | awk '{print $5"\t"$6-1"\t"$7}' \
> $OUTDIR/$PROJECT\_masked_repeats.bed;
