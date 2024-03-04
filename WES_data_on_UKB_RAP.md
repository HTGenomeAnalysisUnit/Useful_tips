# Use WES data on UKB RAP

- [Use WES data on UKB RAP](#use-wes-data-on-ukb-rap)
  - [QC of WES data](#qc-of-wes-data)
    - [QC workflow](#qc-workflow)
    - [Results](#results)
    - [Reference for filtering strategies](#reference-for-filtering-strategies)
    - [Implementation](#implementation)
  - [Perform associations on the WES data](#perform-associations-on-the-wes-data)
    - [Example commands for step2](#example-commands-for-step2)

## QC of WES data

### QC workflow

UKB WES genotypes are provided in RAP as raw calls after merging of sample level data with GLNexus. See the [UKB FAQ](https://www.ukbiobank.ac.uk/media/najcnoaz/access_064-uk-biobank-exome-release-faq_v11-1_final-002.pdf) at point 17 about this.

Based on suggestions from UKB FAQs and internal discussion we decided to apply variant level filters, but not genotypes level filters. The latter would require to re-generate BGEN / VCF file with excessive cost for storage.

1. We split multi-allelic variants in single records and left-align them (bcftools norm). We use the reference genome GRCh38 from [1000G ftp](https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/) as described in the [UKBB documentation](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=914).
2. We compute the following variant level metrics:

| Metric | Description |
|--------|-------------|
| PCT_USABLE_GT | Fraction of GTs that are not missing and with DP >= 10 across all GTs in raw data |
| PCT_GQABOVE20_NOMISS | Fraction of non-missing GTs with GQ >= 20 |

3. We apply a filter on genotypes and we set to missing SNV genotypes with DP < 7 and INDEL genotypes with DP < 10. We then compute the following metrics:

| Metric | Description |
|--------|-------------|
| N_HET | Count of het calls |
| N_BALANCEDHET | Count of het calls with balanced AB (0.2 < AB < 0.8) |
| HWE | HWE test |
| F_MISSING | Fraction of missing calls |
| NS | Number of non-missing calls |

For each WES VCF file in UKB, these metrics are stored in a table with suffix `.metrics.tsv` 

4. Based on these metrics we obtain a list of passing variants that can be used for downstream analysis as follows
   - PCT_USABLE_GT >= 0.9
   - PCT_GQABOVE20_NOMISS >= 0.75
   - N_BALANCEDHET / N_HET >= 0.5
   - F_MISSING < 0.1
   - HWE p-value > 1e-15 (this is not applied for chrX)

The list of passing variant IDs is stored in a file with suffix `.passing_variants.txt`.

### Results

We have applied this filtering strategy to the final release of WES data resulting in the following number of passing variants:

| Chromosome | N_raw_vars | N_pass_vars | PCT_pass_vars |
|------------|------------|-------------|---------------|
| c1         | 2687683    | 2098132     | 0.780647      |
| c2         | 1986463    | 1506609     | 0.758438      |
| c3         | 1572621    | 1201338     | 0.763908      |
| c4         | 1088895    | 799914      | 0.734611      |
| c5         | 1200720    | 909203      | 0.757215      |
| c6         | 1343340    | 1026643     | 0.764247      |
| c7         | 1290956    | 994213      | 0.770137      |
| c8         | 983421     | 752550      | 0.765237      |
| c9         | 1160267    | 904909      | 0.779914      |
| c10        | 1105534    | 837783      | 0.757808      |
| c11        | 1589236    | 1277132     | 0.803614      |
| c12        | 1436008    | 1098130     | 0.76471       |
| c13        | 485362     | 362292      | 0.746437      |
| c14        | 840039     | 654900      | 0.779607      |
| c15        | 936837     | 718983      | 0.767458      |
| c16        | 1300373    | 1059730     | 0.814943      |
| c17        | 1565173    | 1270260     | 0.811578      |
| c18        | 433154     | 326212      | 0.753109      |
| c19        | 1791980    | 1510448     | 0.842893      |
| c20        | 686927     | 552480      | 0.804278      |
| c21        | 289754     | 223251      | 0.770485      |
| c22        | 613857     | 497224      | 0.81          |
| cX         | 652043     | 396406      | 0.607945      |
| total      | 27040643   | 20978742    | 0.775823      |

No QC is provided for variants on chrY due do generally noisy data on chrY.

### Reference for filtering strategies

1. QC from [original UKB WES paper](https://www.nature.com/articles/s41588-021-00885-0#Sec7):

>individual and variant missingness <10%, Hardy Weinberg Equilibrium p-value >10-15, minimum read coverage depth of 7 for SNVs and 10 for indels, at least one sample per site passed the allele balance threshold >0.15 for SNVs and 0.20 for indels.

2. QC from [recent 450k UKB WES](https://www.nature.com/articles/s41586-021-04103-z):

>SNV genotypes with read depth (DP) less than 7 and indel genotypes with read depth less than 10 are changed to no-call genotypes. After the application of the DP genotype filter, a variant-level allele-balance filter is applied, retaining only variants that meet either of the following criteria: (i) at least one homozygous variant carrier; or (ii) at least one heterozygous variant carrier with an allele balance (AB) greater than the cut-off (AB ≥ 0.15 for SNVs and AB ≥ 0.20 for indels).

### Implementation

With this command we split multi-allelic variants in single records and then return a table of metrics for all variants.

```bash
INPUT_PREFIX=$1
REF_GENOME="/mnt/project/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

bcftools norm -m- -f ${REF_GENOME} -Ou ${INPUT_PREFIX}.bcf \
| bcftools +fill-tags -Ou -- -t 'PCT_USABLE_GT=F_PASS(GT!="mis" & FORMAT/DP >= 10)','PCT_GQABOVE20_NOMISS=N_PASS(GT != "mis" & FORMAT/GQ >= 20)/N_PASS(GT != "mis")' \
| bcftools filter -Ou -e '(TYPE = "snp" & FORMAT/DP < 7) | (TYPE = "indel" & FORMAT/DP < 10)' --set-GTs . \
| bcftools +fill-tags -Ou -- -t 'N_HET=N_PASS(GT="het")','N_BALANCEDHET=N_PASS(GT="het" & AD[:1] / DP >= 0.25 & AD[:1] / DP <= 0.75)',F_MISSING,NS,HWE \
| bcftools query -H \
-f '%CHROM:%POS:%REF:%ALT\t%FILTER\t%PCT_USABLE_GT\t%PCT_GQABOVE20_NOMISS\t%N_HET\t%N_BALANCEDHET\t%F_MISSING\t%NS\t%HWE\n' > ${INPUT_PREFIX}.metrics.tsv
```

From the metrics table we can get a list of QCed variants with this command:

```bash
#For autosomes
tail -n+2 ${INPUT_PREFIX}.metrics.tsv \
| awk -F"\t" '{if ($5 == 0) f_het=1; else f_het=$6/$5}; ($2 == "PASS" || $2 == ".") && $3 >= 0.9 && $4 >= 0.75 && $7 < 0.1 && $9 > 1e-15 && f_het >= 0.5' \
| cut -f1 > ${INPUT_PREFIX}.pass_variants.txt

#For X chromosome we do not consider HWE test
tail -n+2 ${INPUT_PREFIX}.metrics.tsv \
| awk -F"\t" '{if ($5 == 0) f_het=1; else f_het=$6/$5}; ($2 == "PASS" || $2 == ".") && $3 >= 0.9 && $4 >= 0.75 && $7 < 0.1 && f_het >= 0.5' \
| cut -f1 > ${INPUT_PREFIX}.pass_variants.txt
```

## Perform associations on the WES data

Here we assume you want to use REGENIE to perform association / rare-variant tests on the WES data.

In step1, you can use the small set of independent variant used normally for GWAS analysis.

When performing step2 one should include the following adjustments:

- exclude chrY from your analysis is reccomended
- use the list of passing variants to automatically subset input variants using the `--extract` option in regenie. The `exome_qc/pass_variants` folder contains a list of pass variants for each chunk. The lists of passing variants for each chromosome are in `exome_qc/pass_variants/per_chromosome`. These files have the same name of the corresponding BGEN/PLINK/VCF file with the suffix `.pass_variants.txt`.
- remove samples with sex discordance using the `--remove` option in regenie. The list of individuals to remove is in the file `exome_qc/ukbb_WES_sex_mismatch_sample_ids.tsv`.
- [UKB best-practises](https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=914) strongly suggests to include a batch covariate when running analysis on WES variants to compensate for different capture kits used across different batches. In particular, samples processed in the first 50k release used a different capture kit than the rest. This information is contained in the [category 170](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=170) of UKB data. The batches for each samples are available in the `exome_qc/ukbb_WES_tranche_release_cov.tsv` file and you should add this to your covariates when running a rare-variant / gene burden test or any other test using WES variant data.

### Example commands for step2

```bash
CHROM=1
INPUT_PREFIX="ukb23157_"
BGEN_FOLDER='/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ BGEN\ format\ -\ final\ release'
QC_FOLDER='/mnt/project/exome_qc'

CHROMOSOME_PREFIX="${INPUT_PREFIX}_c${CHROM}_b0_v1"

regenie \
  --step 2 \
  --bgen ${BGEN_FOLDER}/${CHROMOSOME_PREFIX}.bgen  \
  --sample ${BGEN_FOLDER}/${CHROMOSOME_PREFIX}.sample \
  --phenoFile pheno.tsv \
  --covarFile covar.tsv \
  --extract ${QC_FOLDER}/pass_variants/per_chromosome/${CHROMOSOME_PREFIX}.pass_variants.txt \
  --remove ${QC_FOLDER}/ukbb_WES_sex_mismatch_sample_ids.tsv
  [other options...]
```
