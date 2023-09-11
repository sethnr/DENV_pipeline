#!/usr/bin/env bash

cores=1
while getopts "n:p:s:e:c:d:t:C:O:I:L:" OPTION; do
    case $OPTION in
    n) fname=$OPTARG    ;;
    p) primer_dir=$OPTARG ;;
    
    #python helper scripts
    s) serotype_caller=$OPTARG   ;;
    e) empty_file_maker=$OPTARG   ;;
    
    #mash output
    c) precalls=$OPTARG ;;
    
    #ivar thresholds
    d) depth=$OPTARG ;;
    t) threshold=$OPTARG ;;

    #threading / outdirs
    C) cores=$OPTARG ;;
    O) outdir=$OPTARG ;;
    I) indir=$OPTARG ;;
    L) log=$OPTARG ;;
    *)  echo "option not recognised"
        exit 1
        ;;
    esac
done
shift $(expr $OPTIND - 1 )

read1=$1;
read2=$2;

if [ -z "$precalls" ]; then
    precalls = "${primer_dir}/refs.txt"} 
fi


while IFS= read -r virustype || [[ -n "$virustype" ]]; do 

    fasta=${primer_dir}/${virustype}.fasta
    bed=${primer_dir}/${virustype}.bed
    trimbed=${primer_dir}/${virustype}.trim.bed
    consensus_name=${fname}.${virustype}

    echo "----->>>>>Calling variants against "${virustype}" reference sequence"

    echo "----->>>>>Generating consensus sequence"
    samtools mpileup -aa --reference ${fasta} -A -d 10000 -Q 0 ${indir}/${fname}.${virustype}.sort.bam | ivar consensus -t ${threshold} -m ${depth} -p ${outdir}/${fname}.${virustype}.cons -i ${consensus_name} >> ${log} 2>&1
    
    echo "----->>>>>Aligning consensus cps sequence against the reference serotype "${virustype}" cps sequence"
    nextalign run  --reference ${fasta} --output-fasta ${outdir}/${fname}.${virustype}.out.aln ${outdir}/${fname}.${virustype}.cons.fa >> ${log} 2>&1
    
    if [ -s ${outdir}/${fname%.*}.${virustype}.out.aln ]; then
        echo "----->>>>>Aligning with nextalign successful"
    else
        echo "----->>>>>Aligning with mafft (nextalign not successful)" 
        mafft --quiet --6merpair --keeplength  --addfragments ${outdir}/${fname}.${virustype}.cons.fa ${fasta} > ${outdir}/${fname}.${virustype}.alignment_intermediate.fasta
        grep -A 30000000 ">${consensus_name}" ${outdir}/${fname}.${virustype}.alignment_intermediate.fasta > ${outdir}/${fname}.${virustype}.out.aln
    fi

    echo "----->>>>>Calculating percentage coverage against serotype "${virustype}" cps reference sequence"
    which python
    if [ -s ${trimbed} ]; then
        $CONDA_PREFIX/bin/python ${serotype_caller} --alignment ${outdir}/${fname}.${virustype}.out.aln --bed-file ${trimbed} --outfile ${outdir}/${fname}_all_virustype_info.txt
    else
        $CONDA_PREFIX/bin/python ${serotype_caller}  --alignment ${outdir}/${fname}.${virustype}.out.aln --outfile ${outdir}/${fname}_all_virustype_info.txt
    fi

    echo "----->>>>>Identifying variants"
    samtools mpileup -aa --reference ${fasta} -A -d 0 -Q 0 ${indir}/${fname}.${virustype}.sort.bam | ivar variants -p ${outdir}/${fname}.${virustype}.variants -q 20 -t 0.03 -r ${fasta} >> ${log}     
    echo "----->>>>>>Getting depths"
    bedtools genomecov -d -ibam ${indir}/${fname}.${virustype}.sort.bam > ${outdir}/${fname}.${virustype}.depth.txt; 

    echo "--->>>>> Finished"

done < $precalls
