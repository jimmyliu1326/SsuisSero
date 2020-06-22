#!/usr/bin/env bash

usage(){
  echo "
Usage $0
  Arguments:
  -i|--input    input fastq file
  -o|--outdir   output directory
  -t|--threads  threads
  -h|--help     display help message
"

}

n_threads=4

# parse arguments
if [ $# == 0 ]
then
  usage
  exit 0
fi

opts=`getopt -o hi:o:t: -l input:,outdir:,threads:,help -- "$@"` 
eval set -- "$opts"

while true; do
  case "$1" in
    -i|--input) in_file=$2; shift 2 ;;
    -o|--outdir) out_dir=$2; shift 2 ;;
    -t|--threads) n_threads=$2; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit ;;
  esac
done

# get script dir
Script_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# initialize variables
filename=$(echo $(basename ${in_file%.*}))
ref_path=$(echo $(dirname $Script_DIR)/database/Ssuis_cps2K.fasta)

# main
main() {

  # create outdir
  mkdir -p $out_dir

  echo "Analyzing: $(echo $(basename $in_file))..."

  # run alignment
  minimap2 -ax map-ont -t $n_threads $ref_path $in_file > $out_dir/$filename.sam
  
  # variant calling
  samtools faidx $ref_path
  samtools view -bS -@ $n_threads $out_dir/$filename.sam > $out_dir/$filename.bam
  samtools sort $out_dir/$filename.bam > $out_dir/$filename.sorted.bam
  freebayes -p 1 -f $ref_path $out_dir/$filename.sorted.bam | vcffilter -f "QUAL > 20" > $out_dir/$filename.vcf
  
  #samtools mpileup -B -R -v -A -Q10 -f $ref_path -d 10000 -o $out_dir/$filename.vcf $out_dir/$filename.sorted.bam
  #bcftools view --threads 16 -v snps $out_dir/$filename.vcf > $out_dir/$filename.raw.vcf
  #vcfutils.pl varFilter -Q 20 $out_dir/$filename.raw.vcf > $out_dir/$filename.vcf

  # clean up tmp files
  rm $(echo $(dirname $ref_path)/*.fai)
  rm $out_dir/*.bam
  rm $out_dir/*.sam

}

main
