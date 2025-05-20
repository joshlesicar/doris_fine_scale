# Data processing, read mapping and variant calling  

Initial quality of reads was checked using *FASTQC* v0.12.1 and collated using 
*MULTIQC* v1.18.  

``` bash
module load fastqc
mkdir -p fastqc_output
fastqc *.gz --outdir "./fastqc_output"

module load multiqc
mkdir -p multiqc_output
multiqc . --outdir multiqc_output
```
