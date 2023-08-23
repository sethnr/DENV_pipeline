#!/bin/bash

READS=10000
PROB=1e-100
BLOOM=10
GSIZE=11k

while getopts "p:b:r:" OPTION; do
    case $OPTION in
    p) PROB=$OPTARG    ;;
    r) READS=$OPTARG   ;;
    b) BLOOM=$OPTARG   ;;
    g) GSIZE=$OPTARG   ;;
    T) TEMPDIR=$OPTARG ;;
    *)  echo "option not recognised"
        exit 1
        ;;
    esac
done
shift $(expr $OPTIND - 1 )

let "RLINES = $READS / 2 * 4"

for FQ in $1/*fastq.gz; do
    echo gunzip -dc $FQ \| head -n $RLINES \> ${tempdir}/$(basename ${FQ/fastq.gz/head.fastq})
    gunzip -dc $FQ | head -n $RLINES > ${tempdir}/$(basename ${FQ/fastq.gz/head.fastq})
done

cat ${tempdir}/*head.fastq >  ${tempdir}/${READS}.fastq
rm ${tempdir}/*head.fastq


mash  dist  -m $BLOOM -r -g $GSIZE DENV_all.msh ${1}/head_${READS}.fastq
