# Ssuis_Sero

## Description
This pipeline is designed to rapidly infer Streptococcus suis serotype from Oxford Nanopore Data by first assemblying a draft genome using Miniasm followed by genome polishing with racon and medaka. The processed assembly is subsequently queried against the Cps Blast Database to determine isolate serotype. An additional variant calling step is employed to resolve serotype 2 and 1/2, as well as 1 and 14.

## Usage
<pre>
<b>Required arguments:</b>
-i input raw reads
-o  path to output directory
-s  sample name

<b>Optional arguments:</b>
-h|--help       display help message
-t|--threads    number of threads [4]

_Example Command Line:_
./Ssuis_Sero.sh -s Sample_1 -i /path/to/Sample_1.fastq -o /path/to/output
</pre>

## Conda Installation

### Dependencies
* Miniasm v0.3_r179
* Medaka v1.0.1
* Racon v1.4.13
* Freebayes v1.3.2
* Samtools v1.9
* Minimap2 v2.17
* Vcflib v1.0.0_rc2
