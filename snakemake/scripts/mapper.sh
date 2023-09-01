#!/usr/bin/env bash

# fname=$1
# read1=$2
# read2=$3
# primer_dir=$4
# serotype_caller=$5
# empty_file_maker=$6
# depth=$7
# threshold=$8
# tempdir=$9
# log=$10
cores=1
while getopts "n:p:s:c:e:d:t:C:T:L:" OPTION; do
    case $OPTION in
    n) fname=$OPTARG    ;;
    p) primer_dir=$OPTARG   ;;
    s) serotype_caller=$OPTARG   ;;
    c) precalls=$OPTARG ;;
    e) empty_file_maker=$OPTARG   ;;
    
    C) cores=$OPTARG ;;
    T) tempdir=$OPTARG ;;
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

    if ! [ -f ${primer_dir}/${virustype}.fasta.ann ]; then
        echo "----->>>> making index files for ${virustype}"
        bwa index ${fasta}
    fi

    echo "----->>>>>Mapping reads against serotype "${virustype}" reference sequence"
    bwa mem -v 1 -t ${cores} ${fasta} $read1 $read2 | samtools view -bS -F 4 -F 2048 | samtools sort -@ ${cores} -o ${tempdir}/${fname}.${virustype}.bam >> ${log} 2>&1

    echo "----->>>>>Trimming bam file"
    ivar trim -e -i ${tempdir}/${fname}.${virustype}.bam -b ${bed} -p ${tempdir}/${fname}.${virustype}.trimmed.bam >> ${log} 2>&1

    if ! [ -s  ${tempdir}/${fname%.*}.${virustype}.trimmed.bam ]; then
        echo "no trimmed bam file found, likely because no reads mapped successfully, exiting shell script"
        python ${empty_file_maker} --depth ${depth} --tempdir ${tempdir} --sample-name ${fname} --virus-type ${virustype}
        touch ${tempdir}/${fname}.${virustype}.depth.txt
        continue
    fi

    echo "----->>>>>Sorting bam file"
    samtools sort -@ ${cores} ${tempdir}/${fname}.${virustype}.trimmed.bam -o ${tempdir}/${fname}.${virustype}.sort.bam >> ${log} 2>&1

    echo "----->>>>>Indexing bam file"
    samtools index -@ ${cores} ${tempdir}/${fname}.${virustype}.sort.bam >> ${log} 2>&1

    echo "--->>>>> Mapping Complete"

done < $precalls
