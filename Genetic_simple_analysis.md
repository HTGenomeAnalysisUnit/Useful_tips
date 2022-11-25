# Simple genetic analysis

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
