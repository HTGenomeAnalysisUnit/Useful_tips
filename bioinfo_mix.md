# Bioinfo mix

Random bits of bioinformatics code that I use often

## Estimate kinship from VCF

This uses bcftools and king tools to estimate kinship directly from a VCF/BCF

```bash
#!/bin/bash

inputfile=$1 #Can be VCF or BCF

bcftools norm -Ou -m -any $inputfile \ # Split multi-allelic alleles
    | bcftools norm -Ou -f /resources/hg19/ucsc.hg19.fasta \ # Normalize (optional)
    | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \ # Replace IDs with unique ID
    | plink --bcf /dev/stdin --keep-allele-order --double-id --allow-extra-chr 0 --geno 0.1 --make-bed --out variants # Convert to bed/bim/fam

king -b variants.bed --kinship --prefix kinship  
```

## Convert GTF to refFlat

Picard collect metrics requires a refFlat files of gene annotations. This can be generated from a GTF file using the `gtfToGenePred` tool of UCSC. This can be installed from conda `conda install -c bioconda -c conda-forge ucsc-gtftogenepred` 

Then the following command converts a GTF file to refFlat:

```bash
gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    ${gtf_file} \
    /dev/stdout | \
    awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'
```

## Make introns BED file from GTF

For this we use GTFTOOLS from genemine. This can be downloaded from: https://www.genemine.org/gtftools.php

```bash
gtftools.py -i introns.bed input.gtf
```

## Get contig sizes for GRCh38

GRCh38 assembly contains various categories of contigs. The following command can be used to get the contig sizes, the total size and the size of the effective contigs (excluding decoys and alts) from the fasta file. The HLA contigs are also excluded since they do not represent additional sequences but rather alternative alleles.

```bash
GENOME_NAME="iGenomes_2023.1_GATK_GRCh38"
REF_GENOME="/processing_data/reference_datasets/iGenomes/2023.1/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta" 

grep "^>" $REF_GENOME | tr -s " " "\t" | cut -f1,4 | sed -e 's/>//g' -e 's/LN://g' | grep -v "^HLA" \
| awk '{tot_sum+=$2}; $1 !~ /(_alt|_decoy)/ {sum_effective += $2}; {print ;}; END {print "effective_size", sum_effective; print "total_size", tot_sum}' \
> ${GENOME_NAME}.sizes.txt
```
