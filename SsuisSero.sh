#!/usr/bin/env bash

usage() {
    echo "
Usage: $0

Required arguements:
-i  input file
-o  path to output directory
-s  sample name
-x  input type [fasta or fastq]

Optional arguments:
-h|--help       display help message
-t|--threads    number of threads [Default: 4]
"
}

# default parameters
n_threads=4
sample_name=""
pipeline_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# parse arguments
if [ $# == 0 ]
then
    usage
    exit 0
fi

opts=`getopt -o hi:o:s:t:x: -l help,threads: -- "$@"`
eval set -- "$opts"

while true; do
    case "$1" in
        -i) input_path=$2; shift 2 ;;
        -o) out_dir=$2; shift 2 ;;
        -s) sample_name=$2; shift 2 ;;
        -x) input_type=$2; shift 2 ;;
        -t|--threads) n_threads=$2; shift 2 ;;
        --) shift; break ;;
        -h|--help) usage; exit 0;;
    esac
done

# check if required arguments present
if [ -z $sample_name ]; then
    usage
    echo "Required argument -s missing, exiting"
    exit 1
elif [ -z $out_dir ]; then
    usage
    echo "Required argument -o missing, exiting"
    exit 1
elif [ -z $input_type ]; then
    usage
    echo "Required argument -x missing, exiting"
    exit 1
elif [ -z $input_path ]; then
    usage
    echo "Required argument -i missing, exiting"
    exit 1
fi

# Rapid Assembly

assembly() {

    # set up file structure
    mkdir -p $2
    # Check file integrity

    if test -f $1; then

        # Flye assembly
        flye -t $n_threads --nano-raw $1 -g 2.1m -i 2 -o $2

        # genome polish
        medaka_consensus -t $n_threads -i $1 -d $2/assembly.fasta -o $out_dir -f
        mv $out_dir/consensus.fasta $out_dir/$sample_name.fasta
    else
        echo "$1 cannot be found"
        exit 1
    fi

}

blast_search() {

    if test -f $1; then
        # set up file structure
        mkdir -p $2

        # blast
        blastn -query $1 -db $pipeline_dir/database/Ssuis_Serotyping_blast.db -out $2/blast_res.out -num_threads $n_threads -outfmt 11 -evalue 1.0e-20
    
        # parse blast results
        blast_formatter -archive $2/blast_res.out -outfmt "7 qacc sacc evalue qstart qend sstart send" | awk '!/#/{print}' > $2/blast_res.tab

        if [ $(wc -l < $2/blast_res.tab) -eq 0 ]; then
            serotype="No Hits"
        elif [ $(wc -l < $2/blast_res.tab) -ge 2 ]; then
            serotype=$(cat $2/blast_res.tab | cut -f2 | awk '{gsub(/cps-/,"")}1' | sort -u | paste -sd "|" -)
        else
            serotype=$(cat $2/blast_res.tab | cut -f2 | awk '{gsub(/cps-/,"")}1')
        fi

    else
        echo "$1 cannot be found: Check for errors during assembly step"
        exit 1
    fi    
}

variant_calling() {

    local IFS="|"
    serotypeArray=($serotype)

    if [[ " ${serotypeArray[@]} " =~ " 1 " ]]; then
        mkdir -p $2

        if [[ $input_type == "fastq" ]]; then
            $pipeline_dir/src/ssuis_cpsk_SNP_calling.sh -i $1 -o $2 -t $n_threads
            variants=$(awk '!/#/' $2/*.vcf | cut -f2)
        else
            # genome alignment
            nucmer --maxmatch -b 200 -c 65 -d 0.12 -g 90 -l 20 \
                $pipeline_dir/database/Ssuis_cps2K.fasta $1 -p $2/${sample_name}
            # filter alignments
            delta-filter -1 $2/${sample_name}.delta > $2/${sample_name}.deltafilter
            # call snps
            show-snps -CT $2/${sample_name}.deltafilter | tail -n +5 | cut -f1,2,3 > $2/${sample_name}.vcf
            variants=$(cat $2/*.vcf | cut -f1)
        fi


        # identify variants at position 483
        if echo $variants | grep -q 483; then
            for (( i=0; i<${#serotypeArray[@]}; i++ )); do
                if [[ ${serotypeArray[i]} == "1" ]]; then
                    serotypeArray[i]="14"
                elif [[ ${serotypeArray[i]} == "2" ]]; then
                    serotypeArray[i]="1/2"
                fi
            done
        fi

    elif [[ " ${serotypeArray[@]} " =~ " 2 " ]]; then
        mkdir -p $2
        
        if [[ $input_type == "fastq" ]]; then
            $pipeline_dir/src/ssuis_cpsk_SNP_calling.sh -i $1 -o $2 -t $n_threads
            variants=$(awk '!/#/' $2/*.vcf | cut -f2)
        else
            # genome alignment
            nucmer --maxmatch -b 200 -c 65 -d 0.12 -g 90 -l 20 \
                $pipeline_dir/database/Ssuis_cps2K.fasta $1 -p $2/${sample_name}
            # filter alignments
            delta-filter -1 $2/${sample_name}.delta > $2/${sample_name}.deltafilter
            # call snps
            show-snps -CT $2/${sample_name}.deltafilter | tail -n +5 | cut -f1,2,3 > $2/${sample_name}.vcf
            variants=$(cat $2/*.vcf | cut -f1)
        fi

        # identify variants at position 483
        if echo $variants | grep -q 483; then
            for (( i=0; i<${#serotypeArray[@]}; i++ )); do
                if [[ ${serotypeArray[i]} == "2" ]]; then
                    serotypeArray[i]="1/2"
                elif [[ ${serotypeArray[i]} == "1" ]]; then
                    serotypeArray[i]="14"
                fi
            done
        fi
    
    fi
    
    serotype="$(printf $"|%s" "${serotypeArray[@]}")"
}


write_file() {

    header=$(echo -e "Sample_Name\tSerotype")
    contents=$(echo -e "$sample_name\t${serotype:1}")

    echo $header > $1/${sample_name}_serotyping_res.tsv
    echo $contents >> $1/${sample_name}_serotyping_res.tsv
}

clean() {

    rm $1/*.bam
    rm $1/*.bam.bai
    rm $1/*.hdf
}

# main
main() {

    if [[ $input_type == "fastq" ]]; then
        assembly $input_path $out_dir/assembly
        blast_search $out_dir/${sample_name}.fasta $out_dir/blast_res
        variant_calling $input_path $out_dir/variant_calling
        clean $out_dir
    else
        blast_search $input_path $out_dir/blast_res
        variant_calling $input_path $out_dir/variant_calling
    fi

    write_file $out_dir
    
    echo "Pipeline Finished!"

}

main
