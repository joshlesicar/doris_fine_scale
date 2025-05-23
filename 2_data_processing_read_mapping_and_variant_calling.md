# Data processing, read mapping and variant calling  

## Tools & manuals  
*Fastqc* v0.12.1 [manual.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)  
*MultiQC* v1.18 [manual.](https://github.com/MultiQC/MultiQC)  
*fastp* v0.23.4 [manual.](https://github.com/OpenGene/fastp)  
*kraken2* v2.1.3 [manual.](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)  
*bwa* v0.7.17 [manual.](https://bio-bwa.sourceforge.net/bwa.shtml)  
*samtools* v1.20 [manual.](https://www.htslib.org/doc/samtools.html)  
*picard* v3.1.1 [manual.](https://broadinstitute.github.io/picard/command-line-overview.html#Overview)  
*bcftools* v1.19 [manual.](https://samtools.github.io/bcftools/bcftools.html) 

## Methodology  

### Initial quality check  

Initial quality of reads was checked using *FastQC* v0.12.1 and collated using 
*MultiQC* v1.18. These tools were employed throughout data processing to evaluate the affect of a tool on the reads.  
**Code**  
``` bash
module load fastqc
mkdir -p fastqc_output
fastqc *.gz --outdir "./fastqc_output"

module load multiqc
mkdir -p multiqc_output
multiqc . --outdir multiqc_output
```  
### Adapter contamination removal  
*fastp* v0.23.4 was used to remove adapter contamination, discard reads less than 150 base pairs in length and of Phred score ≤ 20, ensuring a base call precision of 99%.  
**Parameters used**  
```-j```: Output for json format report.  
```-h```: Output for html format report.  
```-q 20```: Discards reads with phred score ≤ 20.  
```-l 150```: Discards reads below 150 bp in length.   
```-detect_adapter_for_pe```: Enables adapter sequence detection for paired ends reads.  
**Code**  
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
### Microorganism and human contamination removal  
*Kraken2* v2.1.3 was then employed to identify potential contaminants from a database of reference sequences including archea, bacteria, viral, plasmid, vector and human. Unclassified reads are isolated for downstream analysis.
Unclassified outputs became ".fq" (uncompressed) files.  
**Parameteres used**  
**Code**  
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
### Read alignment and sorting  
Reads were then aligned to consensus loci ("$CONSENSUS_REF") discovered during initial ddRADseq step (see [ddRADseq and loci discovery](1_ddRADseq_and_loci_discovery.md))
using *bwa* v0.7.17. The maximal exact matches (mem) algorithim was used, as it is recommended for short reads >70 bp. The consenus FASTA was indexed before being used for allignment.To reduce intermediate file creation, reads were also sorted into their
respective genomic coordinates with *samtools* v1.20.  
**Parameters used**  
**Code**  
```bash
cd $CONSENSUS_REF
module load bwa
bwa index dor29_consensus_7858_loci.fasta
``` 
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
### Quality of alignment  
To assess the alignment quality for each sample, samtools function 'flagstat' was employed.  
**code**  
``` bash
module load samtools

for sample in "$INPUT_DIR"/*_sort.bam
do
	samtools flagstat "$sample" > ${sample%_sort.bam}_flagstat.txt
done
```  
### Addition of read groups
Before proceeding with PCR duplication removal, *picard* needs to know read group fields, 
specifically the DNA preparation library identifier (see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)).
As avidity sequencing was used, this field was not included within the fastq metadata field
and must be inserted manually. As only one DNA preperation library was used per sample, 
all fields can use the sample name as a base. This was completed using *picard* v3.1.1 and the 
"AddOrReplaceReadGroups" function.  
**Parameters used**  
**Code**  
```bash
for sample in "$INPUT_DIR"/*.bam
do
base=$(basename "$sample" "_sort.bam")

#Output BAM file with rg tag
OUT="$OUTPUT_DIR/${base}.rg.bam"

echo "Processing sample: $base"

#One library per sample, use sample name as a base for all 5 fields
#avidity sequencing used, but compatible with illumina-style reads
java -jar $PICARD_DIR AddOrReplaceReadGroups \
	I="$sample" \
	O="$OUT" \
	RGID="${base}_id" \
	RGLB="${base}_lib" \
	RGPL=ILLUMINA \
	RGPU="${base}_unit" \
	RGSM="$base"
	
done
```  
### PCR duplicate removal
With the library read group added, PCR duplicate removal using *picard* v3.1.1 can 
proceed. Max file reads was set to slightly below the computing resource maximum of 1024.  
**Parameters used**  
**Code**  
```bash
module load java
module load picard

for sample in $INPUT_DIR/*.rg.bam
do
base=$(basename "$sample" ".rg.bam")

 echo "Processing sample: $base"

java -jar $PICARD_DIR MarkDuplicates \
	I=$INPUT_DIR/${base}.rg.bam \
	O=$OUTPUT_DIR/${base}.bam \
	REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=SILENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=950 \
	M=$REPORT_DIR/${base}_report.txt
done
```  
### Index bam files
Before variant calling, the ```.bam``` files are indexed as *bcftools* (the program used for variant calling) 
requires ```.bai``` index files. This is completed using *samtools* v1.20.  
**Code**  
```bash
module load samtools

for bamfile in $INPUT_DIR/*.bam
do
	echo "Indexing $bamfile"
	samtools index "$bamfile"
done

```  
### Variant calling  
Variant calling was completed with *bcftools* v1.19. Sequences in each ```.bam``` file are 
concatenated together in a mpileup file, which from variants are called using the ```call``` function. 
Both variant and invariant sites are called, as both are required for many population genetic statistics.  
**Parameters used**  
**Code**  
```bash
module load bcftools

bcftools mpileup --threads 15 -a AD,DP,SP -Ou -f $CONSENSUS_REF *.bam | \
bcftools call -f GQ,GP -mv -Oz -o $OUTPUT_DIR/doris29_raw.vcf.gz --threads 15
```

