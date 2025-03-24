#!/bin/bash
# Use the current directory for output
#$ -cwd
#
#
# Request 140G of RAM ... 16*8.75=140
#$ -l h_vmem=30G
#
# Request the large mem queue
#$ -q 4-day
#$ -j y
# Request stack hard limit (for OpenMP)
#$ -l h_stack=5M
#
# Request local tmp space
##$ -l loc2tmp=40G
#
# Others
#$ -notify

##
# Working skeleton for clean-up pipeline

set -e



#########################################################################
# Parameters to set:
#
# MYDIR =  directory to run the analysis
# SAMPLEFILE = tsv of input filenames. 
#              Col1: sample ID
#              Col2: R1 fastq filename 
#              Col3: R2 fastq filename
# OUTPUT = name of output directory
# PROTEIN = path to humann3 protein database (e.g., /data/uniref)
# CHOCO = path to humann3 chocophlan database (e,g., /data/chocophlan)
# HUMAN = path to human genome bloom filter
# MOUSE = path to mouse genome bloom filter
# ALL_RRNA = path to sortmerna reference fasta
# K2DB = path to kraken2 database
#########################################################################

MYDIR=processing/
SAMPLEFILE=$MYDIR/files.txt
SAMPLEDIR=raw_fastq/
OUTPUT=output/
PROTEIN=db/humann/data/uniref
CHOCO=db/humann/data/chocophlan
HUMAN=db/biobloom23/human38.bf
MOUSE=db/biobloom23/GRCm49.bf
ALL_RRNA=db/rsortmerna4/all-rRNA.fasta
K2DB=db/kraken2/gtdb



SAMPLE=$(sed -n ${SGE_TASK_ID}p $SAMPLEFILE | cut -f1)
SAMPLE_R1=$(sed -n ${SGE_TASK_ID}p $SAMPLEFILE | cut -f2)
SAMPLE_R2=$(sed -n ${SGE_TASK_ID}p $SAMPLEFILE | cut -f3)

#die function
warn () {
    echo "$0:" "$@" >&2
}
die () {
    rc=$1
    shift
    warn "$@"
    exit $rc
}

isitthere () {
  [ -s $1 ] || die 1 "File $1 not found or zero sized."
}

#Make the workdir for RCF
#Traps for errors

[[ -d $MYDIR/tmp ]] || mkdir -p $MYDIR/tmp

WORKDIR=$(mktemp -d $MYDIR/tmp/${USER}_tmp.XXXXXXXX)
echo "Workdir: $WORKDIR"

cleanup() {
        #cleanup code for tmp
        /bin/rm -rf $WORKDIR

        exit
}

trap cleanup USR1 USR2 EXIT


export PATH=$HOME/bin:$PATH



cd $WORKDIR

module load bbmap

bbduk.sh threads=8 ref=adapters unpigz=t ktrim=r k=23 mink=6 hdist=1 tpe tbo \
        in=$SAMPLEDIR/$SAMPLE_R1 in2=$SAMPLEDIR/$SAMPLE_R2 interleaved=t \
        out=stdout.fq -Xmx4g |
bbduk.sh threads=8 in=stdin.fq interleaved=t out=stdout.fq qtrim=r trimq=15 -Xmx4g |
bbduk.sh threads=8 k=23 ref=phix hdist=2 minkmerhits=2 unpigz=t \
        in=stdin.fq interleaved=t out=stdout.fq -Xmx4g |
bbduk.sh threads=8 entropy=0.7 entropymask=f pigz=t \
        in=stdin.fq interleaved=t \
        out=${SAMPLE}_R1.bbduk.fastq.gz out2=${SAMPLE}_R2.bbduk.fastq.gz\
        -Xmx4g

module load biobloomtools

biobloomcategorizer -t $NSLOTS -p ${SAMPLE} -f "$HUMAN $MOUSE" \
	-e -i --fq ${SAMPLE}_R1.bbduk.fastq.gz ${SAMPLE}_R2.bbduk.fastq.gz

module load sortmerna

sortmerna --threads $NSLOTS \
	--ref $ALL_RRNA \
	--reads ${SAMPLE}_noMatch_1.fq \
	--reads ${SAMPLE}_noMatch_1.fq \
	--workdir $MYDIR/sortmerna/${SAMPLE}_workdir \ 
	--idx $MYDIR/sortmerna/${SAMPLE}_idx \
	--best 1 \
	--paired_in \
	--fastx \
	--out2 \
	--other \
	$MYDIR/sortmerna/${SAMPLE}.sortmerna.out

reformat.sh in1=$MYDIR/sortmerna/${SAMPLE}.sortmerna.out_fwd.fq in2=$MYDIR/sortmerna/${SAMPLE}.sortmerna.out_rev.fq out=${SAMPLE}.cleaned.fastq.gz
conda deactivate

rm -rf $MYDIR/sortmerna/${SAMPLE}_workdir

export PATH=$HOME/bin:$PATH


module load humann

[[ -d $MYDIR/$OUTDIR ]] || mkdir -p $MYDIR/$OUTDIR

humann3 -vvv --input ${SAMPLE}.cleaned.fastq.gz --input-format fastq.gz \
	--output $MYDIR/$OUTDIR --output-basename ${SAMPLE} \
	--threads $NSLOTS \
	--protein-database $PROTEIN \
	--nucleotide-database $CHOCO

rm ${SAMPLE}_humann_temp/*.sam


module load kraken2/2.0.9 bracken

kraken2 --threads $NSLOTS \
	--db $K2DB \
	--confidence 0.25 \
	--report $WORKDIR/${SAMPLE}.k2.025.report \
	--output ${SAMPLE}.kraken.out \
	--paired $MYDIR/sortmerna/${SAMPLE}.sortmerna.out_fwd.fq $MYDIR/sortmerna/${SAMPLE}.sortmerna.out_rev.fq


bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.phylum.025.bracken  -l P -t 10
bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.class.025.bracken  -l C -t 10
bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.order.025.bracken  -l O -t 10
bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.family.025.bracken  -l F -t 10
bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.genus.025.bracken  -l G -t 10
bracken -d $K2DB -i ${SAMPLE}.k2.025.report -o ${SAMPLE}.species.025.bracken -l S -t 10


[[ -d $MYDIR/kraken_out_cleaned_final ]] || mkdir -p $MYDIR/kraken_out_cleaned_final
mv $WORKDIR/${SAMPLE}.kraken.out $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.k2.025.report $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.phylum.025.bracken $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.class.025.bracken $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.order.025.bracken $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.family.025.bracken $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.genus.025.bracken $MYDIR/kraken_out_cleaned_final/
mv $WORKDIR/${SAMPLE}.species.025.bracken $MYDIR/kraken_out_cleaned_final/
 
echo "Done!"
