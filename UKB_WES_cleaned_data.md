# Clean WES data on UKB RAP

UKB WES genotypes are provided in RAP as raw calls after merging of sample level data with GLNexus. See the [UKB FAQ](https://www.ukbiobank.ac.uk/media/najcnoaz/access_064-uk-biobank-exome-release-faq_v11-1_final-002.pdf) at point 17 about this.

Based on suggestions from UKB FAQs and internal discussion we decided to apply variant level filters, but not genotypes level filters. The latter would require to re-generate BGEN / VCF file with excessive cost for storage.

1. We extract bi-allelic variants only
2. We set GT to missing for variants with DP < 10
3. We compute the following variant level metrics:

| Metric | Description |
|--------|-------------|
| PCTDPABOVE10_ALL | Fraction of GTs with DP >= 10 |
| PCTGQABOVE20_NOMISS | Fraction of non-missing GTs with GQ >= 20 |
| N_HET | Count of het calls |
| NBALANCEDHET | Count of het calls with balanced AB (0.2 < AB < 0.8) |
| HWE | HWE test |
| F_MISSING | Fraction of missing calls |
| NS | Number of non-missing calls |

For each WES VCF file in UKB, these metrics are stored in a table with suffix `.metrics.tsv` and a list of passing variants is stored in a file with suffix `.passing_variants.txt`.

## Reference for filtering strategies

1. QC from [original UKB WES paper](https://www.nature.com/articles/s41588-021-00885-0#Sec7):

>individual and variant missingness <10%, Hardy Weinberg Equilibrium p-value >10-15, minimum read coverage depth of 7 for SNVs and 10 for indels, at least one sample per site passed the allele balance threshold >0.15 for SNVs and 0.20 for indels.

2. QC from [recent 450k UKB WES](https://www.nature.com/articles/s41586-021-04103-z):

>SNV genotypes with read depth (DP) less than 7 and indel genotypes with read depth less than 10 are changed to no-call genotypes. After the application of the DP genotype filter, a variant-level allele-balance filter is applied, retaining only variants that meet either of the following criteria: (i) at least one homozygous variant carrier; or (ii) at least one heterozygous variant carrier with an allele balance (AB) greater than the cut-off (AB ≥ 0.15 for SNVs and AB ≥ 0.20 for indels).

'F_PASS(DP>=10 & GT!="mis")> 0.9'

## Implementation

With this command we return a table of metrics for all bi-allelic variants.

```bash
INPUT_PREFIX=$1

bcftools view -m2 -M2 -Ou ${INPUT_PREFIX}.vcf.gz \
| bcftools filter -Ou -e 'FORMAT/DP < 10' --set-GTs . \
| bcftools +fill-tags -Ou -- -t 'PCTDPABOVE10_ALL=F_PASS(FORMAT/DP >= 10)','PCTGQABOVE20_NOMISS=F_PASS(GT != "mis" & FORMAT/GQ >= 20)','N_HET=N_PASS(GT="het")','NBALANCEDHET=N_PASS(GT="het" & AD[:1] / DP >= 0.25 & AD[:1] / DP <= 0.75)',F_MISSING,NS,HWE \
| bcftools query -H \
-f '%CHROM:%POS:%REF:%ALT\t%FILTER\t%PCTDPABOVE10_ALL\t%PCTGQABOVE20_NOMISS\t%N_HET\t%NBALANCEDHET\t%F_MISSING\t%NS\t%HWE\n' > ${INPUT_PREFIX}.metrics.tsv
```

From the metrics table we can get a list of QCed variants with this command:

```bash
awk -F'\t' '$3 > 0.9 && $4 > 0.9 && $5 > 0 && $6 > 0 && $7 < 0.1 && $8 > 0' ${INPUT_PREFIX}.metrics.tsv | cut -f1 > ${INPUT_PREFIX}.passing_variants.txt
```

## Use with regenie

One can use the list of passing variants to automatically subset input variants in regenie with the `--extract` option at step2 stage. In step1, you can use the small set of independent variant used normally for GWAS analysis.

We have also compiled a list of samples to remove in the file ##### and you should also filter your samples accordingly using the `--remove` option in regenie.

Finally, [UKB best-practises](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=914) strongly suggest to include a batch covariate when running analysis on WES variants to compensate for different capture kits used across different batches.

The batches for each samples are available in this file and should add this to your covariates when running a rare-variant / gene burden test or any other test using WES variant data.
