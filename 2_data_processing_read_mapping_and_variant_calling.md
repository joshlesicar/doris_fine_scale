# Data processing, read mapping and variant calling  

Initial quality of reads was checked using *FastQC* v0.12.1 and collated using 
*MultiQC* v1.18.  

``` bash
module load fastqc
mkdir -p fastqc_output
fastqc *.gz --outdir "./fastqc_output"

module load multiqc
mkdir -p multiqc_output
multiqc . --outdir multiqc_output
```  
*fastp* v0.23.4 was used to remove adapter contamination, discard reads less than 150 base pairs in length and of Phred score $le; 20 (Q20), ensuring a base call precision of 99%.

``` bash
module load fastp

for f1 in *_R1.fastq.gz
do
        f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
		fastp -i $f1 -I $f2 -o "$OUT_DIR/step1.0_${f1}" -O "$OUT_DIR/step1.0_${f2}" \
              -j "$REPORT_DIR/step1.0_${f1}.json" \
              -h "$REPORT_DIR/step1.0_${f1}.html" \
			        -q 20 \
			        -l 150 \
              --detect_adapter_for_pe
done
```