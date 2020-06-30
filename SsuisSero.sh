#!/usr/bin/env bash

usage() {
    echo "
Usage: $0

Required arguements:
-i input raw reads
-o  path to output directory
-s  sample name

Optional arguments:
-h|--help       display help message
-t|--threads    number of threads [4]
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

opts=`getopt -o hi:o:s:t: -l help,threads: -- "$@"`
eval set -- "$opts"

while true; do
    case "$1" in
        -i) read_path=$2; shift 2 ;;
        -o) out_dir=$2; shift 2 ;;
        -s) sample_name=$2; shift 2 ;;
        -t|--threads) n_threads=$2; shift 2 ;;
        --) shift; break ;;
        -h|--help) usage; exit 0;;
    esac
done

if [ -z $sample_name ]; then
    echo "Sample name was not given exiting"
    exit 1
fi

# Rapid Assembly

assembly() {

    # set up file structure
    mkdir -p $2
    # Check file integrity

    if test -f $1; then

        # read overlap
        minimap2 -t $n_threads -x ava-ont $1 $1 > $2/overlaps.paf

        # generate OLC graph
        miniasm -f $1 $2/overlaps.paf > $2/graph.gfa

        # minipolish
        minipolish -t $n_threads --rounds 2 $1 $2/graph.gfa > $2/graph_polished.gfa

        # convert to fasta format
        awk '$1 ~/S/ {print ">"$2"\n"$3}' $2/graph_polished.gfa > $2/graph_polished.fasta
        
        # genome polish
        medaka_consensus -t $n_threads -i $1 -d $2/graph_polished.fasta -o $out_dir -f
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
            serotype=$(cat $2/blast_res.tab | cut -f2 | awk '{gsub(/cps-/,"")}1' |  paste -sd "|" -)
        else
            serotype=$(cat $2/blast_res.tab | cut -f2 | awk '{gsub(/cps-/,"")}1')
        fi

    else
        echo "$1 cannot be found: Check for errors during assembly step"
        exit 1
    fi    
}

variant_calling() {

    if [ $(echo $serotype | grep "2\|1") ]; then

        mkdir -p $2
        $pipeline_dir/src/ssuis_cpsk_SNP_calling.sh -i $1 -o $2 -t $n_threads

        # identify variants at position 483
        variants=$(awk '!/#/' $2/*.vcf | cut -f2)
        if echo $variants | grep -q 483; then
            serotype=$(echo $serotype | awk '{gsub(/1/,"14")}1')
            serotype=$(echo $serotype | awk '{gsub(/2/,"1/2")}1')
            
        fi
    
    fi
            
}

write_file() {

    header=$(echo -e "Sample_Name\tSerotype")
    contents=$(echo -e "$sample_name\t$serotype")

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

    assembly $read_path $out_dir/assembly
    blast_search $out_dir/$sample_name.fasta $out_dir/blast_res
    variant_calling $read_path $out_dir/variant_calling
    write_file $out_dir
    clean $out_dir
    echo "Pipeline Finished!"

}

main
