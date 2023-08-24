#!/bin/bash

READS=10000
PROB=1e-100
BLOOM=10
GSIZE=11k

while getopts "p:b:r:g:T:" OPTION; do
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


TEMPDIR=`mktemp -d`;
sleep 1 

if [[ ! "$TEMPDIR" || ! -d "$TEMPDIR" ]]; then
    echo "!!!!!"
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

echo "mashing ${READS} reads against hashes"
mash  dist  -m $BLOOM -r -g $GSIZE DENV_all.msh ${TEMPDIR}/${READS}.fastq

#rm -r ${TEMPDIR}
