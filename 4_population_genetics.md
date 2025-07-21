
# Population genetics

## Tools & manuals

*vcftools* v0.1.16 [manual.](https://vcftools.github.io/man_latest.html)

## Structure

### Create structure data file using plink

    plink --vcf dataset1.recode.vcf --allow-extra-chr --recode --out dataset1temp
    plink --file dataset1temp --allow-extra-chr --recode structure --out dataset1_structure

### Create parameter files

### Use structure threader to qsub jobs

## Population genetics of largest cluster

### Create new dataset

Dataset 2 will be comprised of all individuals within the largest
cluster. However, removal of individuals can affect filtering steps
(e.g. maf), so filtering steps from and including maximum observed
heterozygostiy will be reapplied following removal of the following
individuals.

    RN02_L0056
    RN02_L0062
    RN02_L0125
    RN02_L0127
    RN02_L0139
    RN02_L0141
    RN02_L0142
    RN02_L0145
    RN02_L0146
    RN02_L0148
    RN02_L0149
    RN02_L0150
    RN02_L0152
    RN02_L0154
    RN02_L0156
    RN02_L0157
    RN02_L0158
    RN02_L0160
    RN02_L0163
    RN02_L0164

A text file, `indv_remove_dataset2` was created with these individuals.
They will be removed using *vcftools* function `--remove`.

    vcftools --vcf max_missing.recode.vcf --remove indv_remove_dataset2.txt --recode --recode-INFO-all --out dataset2_maxmissing

Data loss: After filtering, kept 144 out of 164 Individuals. After
filtering, kept 68504 out of a possible 68504 Sites.

**Filter by maximum observed heterozygosity, minor allele frequency and
keep 1 SNP per locus.** For methodology please refer to [Variant
filtering.](3_variant_filtering.md)

Data loss from heterozygosity: After filtering, kept 60681 out of a
possible 68504 Sites.

Data loss from minor allele frequency: After filtering, kept 11243 out
of a possible 60681 Sites.

Data loss from retaining 1 SNP per locus: After filtering, kept 4218 out
of a possible 11243 Sites.
