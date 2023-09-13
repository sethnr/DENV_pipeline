#!/usr/bin/env bash

cores=1


genome_size="11k"

while getopts "i:o:g:m:C:" OPTION; do
    case $OPTION in
    i) infile=$OPTARG      ;;
    o) outfile=$OPTARG     ;;
    m) mashout=TRUE         ;;
    g) genome_size=$OPTARG ;;    
    C) cores=$OPTARG       ;;
    *)  echo "option not recognised"
        exit 1
        ;;
    esac
done
shift $(expr $OPTIND - 1 )

if [ -z "$mashout" ]; then
    FASTA=$1
    echo "generating bwa index for ${FASTA}"  2>&1
    echo bwa index ${FASTA} 2>&1
    bwa index ${FASTA}
else 
    FASTAS=$@;
    echo "making mash indices for ${FASTAS}"  2>&1
    for $F in ${FASTAS}; do 
        echo mash sketch -g ${genome_size} ${F} 2>&1 
        mash sketch -g ${genome_size} ${F} ; 
    done

    mashidx=("${FASTAS[@]}")

    for i in "${!mashidx[@]}"; do
        echo mashidx[$i]  --> ${mashidx[$i]/.fasta/.msh}
        mashidx[$i]=${mashidx[$i]/.fasta/.msh}
    done

    echo mash paste $mashout ${mashidx}
    mash paste $mashout ${mashidx}
fi