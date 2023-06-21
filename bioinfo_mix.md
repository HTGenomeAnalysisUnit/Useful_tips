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