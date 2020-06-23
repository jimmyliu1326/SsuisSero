# Ssuis_Sero

## Description
This pipeline is designed to rapidly infer Streptococcus suis serotype from Oxford Nanopore Data by first assemblying a draft genome using Miniasm followed by genome polishing with racon and medaka. The processed assembly is subsequently queried against the Cps Blast Database to determine isolate serotype. An additional variant calling step is employed to resolve serotype 2 and 1/2, as well as 1 and 14.

## Usage
```
Required arguements:
-i input raw reads
-o  path to output directory
-s  sample name
Optional arguments:
-h|--help       display help message
-t|--threads    number of threads [4]

Example Command Line:
./Ssuis_Sero.sh -s Sample_1 -i /path/to/Sample_1.fastq -o /path/to/output
``` 

## Conda Installation
