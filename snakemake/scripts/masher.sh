#!/bin/bash

READS=10000
PROB=1e-10
DIST=0.5
BLOOM=10
GSIZE=11k
PREFIX="outfile"
MASHREF="DENV_all.msh"

while getopts "f:p:d:b:r:g:T:o:" OPTION; do
    case $OPTION in
    f) MASHREF=$OPTARG ;;
    p) PROB=$OPTARG    ;;
    d) DIST=$OPTARG    ;;
    r) READS=$OPTARG   ;;
    b) BLOOM=$OPTARG   ;;
    g) GSIZE=$OPTARG   ;;
    o) PREFIX=$OPTARG  ;;
    T) TEMPDIR=$OPTARG ;;
    *)  echo "option not recognised"
        exit 1
        ;;
    esac
done
shift $(expr $OPTIND - 1 )


TEMPDIR=`mktemp -d`;
sleep 1 

if [[ ! "$TEMPDIR" || ! -d "$TEMPDIR" ]]; then
    echo $TEMPDIR not found
    exit 1
fi

echo "getting top ${READS} reads from infiles" 1>&2
let "RLINES = $READS / 2 * 4"
for FQ in $1/*fastq.gz; do
    BASENAME=$(basename ${FQ/fastq.gz/head.fastq})
    echo $BASENAME 1>&2
    echo gunzip -dc $FQ \| head -n $RLINES \> ${TEMPDIR}/${BASENAME}
    gunzip -dc $FQ | head -n $RLINES > ${TEMPDIR}/${BASENAME}
done

cat ${TEMPDIR}/*head.fastq >  ${TEMPDIR}/${READS}.fastq
rm ${TEMPDIR}/*head.fastq

echo "mashing ${READS} reads against hashes" >&2
mash  dist  -m $BLOOM -r -g $GSIZE ${MASHREF} ${TEMPDIR}/${READS}.fastq > ${PREFIX}_mash.txt


#filter for max prob / dist, print genome names
echo "pulling matches below ${DIST} / ${PROB}" >&2
awk -v dist=$DIST -v prob=$PROB '($3+0 < dist+0 && $4+0 < prob+0) {sub(".fasta","",$1); print $1}' ${PREFIX}_mash.txt > ${PREFIX}_calls.txt


echo "rm ${TEMPDIR}" >&2
rm -r ${TEMPDIR}

echo "fin" >&2
