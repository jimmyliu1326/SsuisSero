# SsuisSero

## Description
This pipeline is designed to rapidly infer Streptococcus suis serotype from Oxford Nanopore Data by first assemblying a draft genome using Miniasm followed by genome polishing with racon and medaka. The processed assembly is subsequently queried against the Cps Blast Database to determine isolate serotype. An additional variant calling step is employed to resolve serotype 2 and 1/2, as well as 1 and 14.

## Usage
<pre>
<b>Required arguments:</b>
-i  input raw reads
-o  path to output directory
-s  sample name

<b>Optional arguments:</b>
-h|--help       display help message
-t|--threads    number of threads [4]

<b>Example Command Line:</b>
SsuisSero.sh -s Sample_1 -i /path/to/Sample_1.fastq -o /path/to/output
</pre>

## Conda Installation
The recommended method of installation for this pipeline is with Conda. Set up a new conda environment and run:

```
conda install -c bioconda ssuissero
```

### Dependencies
* Medaka >= 1.0.1
* Freebayes >= 1.3.2
* Samtools >= 1.9
* Vcflib >= 1.0.0_rc2
* Blast >= 2.6
* Mummer = 3.23
