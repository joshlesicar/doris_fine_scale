# Data processing, read mapping and variant calling  

Initial quality of reads was checked using *FastQC* v0.12.1 and collated using 
*MultiQC* v1.18. These tools were employed throughout data processing to evaluate the affect of a tool on the reads.

``` bash
module load fastqc
mkdir -p fastqc_output
fastqc *.gz --outdir "./fastqc_output"

module load multiqc
mkdir -p multiqc_output
multiqc . --outdir multiqc_output
```  
*fastp* v0.23.4 was used to remove adapter contamination, discard reads less than 150 base pairs in length and of Phred score ≤ 20 (Q20), ensuring a base call precision of 99%.

``` bash
module load fastp

for f1 in *_R1.fastq.gz
do
    f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
		fastp -i $f1 -I $f2 -o "$OUT_DIR/${f1}" -O "$OUT_DIR/${f2}" \
    -j "$REPORT_DIR/step1.0_${f1%%_R1.fastq.gz}.json" \
    -h "$REPORT_DIR/step1.0_${f1%%_R1.fastq.gz}.html" \
		-q 20 \
		-l 150 \
    --detect_adapter_for_pe
done
```

*Kraken2* v2.1.3 was then employed to identify potential contaminants from a database of reference sequences including archea, bacteria, viral, plasmid, vector and human. Unclassified reads are isolated for downstream analysis.
Unclassified outputs became ".fq" (uncompressed) files.
``` bash
for R1 in "$INPUT_DIR"/*_R1.fastq.gz
do
    R2=${R1%%_R1.fastq.gz}"_R2.fastq.gz"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    if [[ -f "$R2" ]]; then
        echo "Processing sample: $SAMPLE"
   		kraken2 \
  			--use-names \
			  --threads 12 \
  			--db "$KRAKEN_DATABASE" \
  			--gzip-compressed \
  			--paired "$R1" "$R2" \
  			--classified-out "${CLASSIFIED_DIR}/${SAMPLE}_classified#.fq" \
  			--unclassified-out "${UNCLASSIFIED_DIR}/${SAMPLE}_unclassified#.fq" \
  			--report "${REPORT_DIR}/${SAMPLE}_report.txt"
    fi
done

```  
Reads were then aligned to consensus loci ("$CONSENSUS_REF") discovered during initial ddRADseq step (see [ddRADseq and loci discovery](1_ddRADseq_and_loci_discovery.md))
using *bwa* v0.7.17. To reduce intermediate file creation, reads were also sorted into their
respective genomic coordinates with *samtools* v1.20.
``` bash
module load bwa
module load samtools

# Get individual names by stripping _1.fq
INDS=($(for i in "$SEQUENCES"/*_1.fq; do basename "$i" | sed -E 's/_unclassified.*//'; done))

for IND in "${INDS[@]}"; do
    FORWARD="$SEQUENCES/${IND}_unclassified_1.fq"
    REVERSE="$SEQUENCES/${IND}_unclassified_2.fq"
    OUTPUT="$OUTDIR/${IND}_sort.bam"

    echo "Aligning $IND with bwa"
    bwa mem -M -t 15 "$CONSENSUS_REF" "$FORWARD" "$REVERSE" \
        | samtools view -b -@ 15 \
        | samtools sort -@ 15 -T "$OUTDIR/${IND}" > "$OUTPUT"
done
```