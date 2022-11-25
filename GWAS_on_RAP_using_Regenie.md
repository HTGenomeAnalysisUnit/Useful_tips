# Tutorial for running a GWAS on the RAP using Regenie

- [Tutorial for running a GWAS on the RAP using Regenie](#tutorial-for-running-a-gwas-on-the-rap-using-regenie)
  - [Prerequisite](#prerequisite)
  - [1.Create target cohort](#1create-target-cohort)
  - [2.Prepare input files](#2prepare-input-files)
    - [Set-up your folders](#set-up-your-folders)
    - [Create command to execute](#create-command-to-execute)
    - [Run command in dxpy](#run-command-in-dxpy)
  - [3.Plink QC](#3plink-qc)
  - [4.Regenie step 1](#4regenie-step-1)
  - [5.Regenie step 2](#5regenie-step-2)

## Prerequisite

login in dxpy (DNA Nexus CLI tool used to interact with the platform) by typing

```bash
conda activate dxpy_0.327
dx login
```

## 1.Create target cohort

Since this is only an example, we extracted a subset from the available sample. We restricted the analysis to white British women, between 55 and 60 years old, with not missing BMI values. This subset cohort was created using the [Cohort browser](https://documentation.dnanexus.com/user/cohort-browser#opening-datasets-using-the-cohort-browser) and then exported into text file by using the [Table exporter app](https://ukbiobank.dnanexus.com/app/table-exporter).
Find [here](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-phenotypic-data-as-a-file#selecting-fields-of-interest-in-the-cohort-browser) further details on how to access phenotypic data as a file check.

## 2.Prepare input files

Input file needed:

- Phenotype file reporting, for each individual, values for traits of interest and covariates. Individuals iid should be reported in the first two columns and named "FID" and "IID"
- Genetic data file format. To speed-up the process, here we use genotype calls in Plink format (which need to be provided in a single, genome-wide file)

### Set-up your folders

```bash
dx mkdir -p /Regenie_test/step1/
dx mkdir -p /Regenie_test/step2/
```

### Create command to execute

We are first creating a file listing all plink files (split by chromosome) to merge. Symbolic links are used to bypass the original path where the file are stored, which contains space characters and causes plink to fail (symbolic links are located in folder simply to make things more tidy). Then we create a file listing all samples id to extract from our newly merged plink files - this last step is not necessary for the analysis sake, but only to generate smaller output files

```bash
run_plink_merge="mkdir temp; 
  cd temp;
  ln -s /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* ./;
  cd ..;
  ls temp/*.bed | sed 's/.bed//g' > plink_files_to_merge.txt;
  awk 'NR==1 {print \"FID\",\"IID\",\$2,\$3}; NR > 1 {print \$1, \$1, \$2, \$3}' random_cohort.tsv > regenie_pheno_form_input.tsv
  awk '{print \$1,\$1}' regenie_pheno_form_input.tsv > iid_list_plink.tsv
  plink --merge-list plink_files_to_merge.txt \
    --keep iid_list_plink.tsv \
    --make-bed --autosome --out ukb22418_c1_22_v2_merged_subset
  rm plink_files_to_merge.txt"
```

### Run command in dxpy

We use the Swiss army knife app to run our command, since it includes Plink and Regenie softwares

```bash
dx run swiss-army-knife \
  -iin="/Regenie_test/data/random_cohort.tsv" \
  -icmd="${run_plink_merge}" \
  --tag="pheno_input_and_plink_merge" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="/Regenie_test/data/" --brief --yes
```

## 3.Plink QC

Perform standard QC of newly merged, Plink genotype file

```bash
run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged_subset \
 --autosome --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
 --write-snplist --write-samples --no-id-header --out geno_array_snps_qc_pass"

dx run swiss-army-knife \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.bed" \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.bim" \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.fam" \
 -icmd="${run_plink_qc}" \
 --tag="plink_qc" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Regenie_test/data/" --brief --yes
```

## 4.Regenie step 1

```bash
run_regenie_step1="regenie --step 1 \
 --bed ukb22418_c1_22_v2_merged_subset \
 --phenoFile /mnt/project/Regenie_tes/data/regenie_pheno_form_input.tsv \
 --covarFile /mnt/project/Regenie_tes/data/regenie_pheno_form_input.tsv \
 --extract /mnt/project/Regenie_tes/data/geno_array_snps_qc_pass.snplist \
 --phenoCol p21001_i0 \
 --covarCol p21022 \
 --out bmi_gwas_test \
 --bsize 1000 --lowmem --qt --loocv --threads 16 --gz"

dx run swiss-army-knife \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.bed" \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.bim" \
 -iin="/Regenie_test/data/ukb22418_c1_22_v2_merged_subset.fam" \
 -icmd="${run_regenie_step1}" \
 --tag="regenie_step1" \
 --instance-type "mem1_ssd1_v2_x16" \
 --destination="/Regenie_test/step1/" --brief --yes
```

## 5.Regenie step 2

```bash
for chr in {1..22}
do

  run_regenie_step2="regenie --step 2 \
    --bgen /mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c${chr}_b0_v3.bgen \
    --sample /mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c${chr}_b0_v3.sample \
    --phenoFile /mnt/project/Regenie_test/data/regenie_pheno_form_input.tsv \
    --covarFile /mnt/project/Regenie_test/data/regenie_pheno_form_input.tsv \
    --pred bmi_gwas_test_pred.list \
    --phenoCol p21001_i0 \
    --covarCol p21022 \
    --minMAC 100 \
    --out bmi_gwas_test_chr${chr} \
    --bsize 200 --qt --threads 16 --gz"

    dx run swiss-army-knife \
      -iin="/Regenie_test/step1/bmi_gwas_test_pred.list" \
      -iin="/Regenie_test/step1/bmi_gwas_test_1.loco.gz" \
      -icmd="${run_regenie_step2}" \
      --name="regenie_step2_chr${chr}" \
      --tag="regenie_step2" \
      --instance-type="mem1_ssd1_v2_x16" \
      --destination="/Regenie_test/step2/" --brief --yes
done
```
