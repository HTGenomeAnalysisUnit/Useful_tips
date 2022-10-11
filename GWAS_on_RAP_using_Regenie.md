# Tutorial for running a GWAS on the RAP using Regenie

- [Create target cohort](#create-target-cohort)
- [Prepare_input_files](#prepare-input-files)
- [Plink_QC](#plink-qc)
- [GWAS](#GWAS)
  - [Regenie_step1](#regenie-step1)
  - [Regenie_step2](#regenie-step2)


### 0) Prerequisite - login in dxpy (DNA Nexus CLI tool used to interact with the platform) by typing:
```bash
conda activate dxpy_0.327
dx login
```


### 1) Create target cohort

Since this is only an example, we extracted a subset from the available sample. We restricted the analysis to women, between 55 and 60 years old, with not missing BMI values. This subset cohort was created using the [Cohort browser](https://documentation.dnanexus.com/user/cohort-browser#opening-datasets-using-the-cohort-browser) and then exported into text file by using the [Table exporter app](https://ukbiobank.dnanexus.com/app/table-exporter).
Find [here](https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-phenotypic-data-as-a-file#selecting-fields-of-interest-in-the-cohort-browser) further details on how to access phenotypic data as a file check.



### 2) Prepare input files
Input file needed:
  - Phenotype file reporting, for each individual, values for traits of interest and covariates. Individuals iid should be reported in the first two columns and named "FID" and "IID"
  - Genetic data file format. To speed-up the process, here we use genotype calls in Plink format (which need to be provided in a single, genome-wide file)

#### Set-up your folders:
```bash
dx mkdir -p /Regenie_test/step1/
dx mkdir -p /Regenie_test/step2/
dir="/Regenie_test/data/"
```

#### Create command to execute:
We are first creating a file listing all plink files (split by chromosome) to merge. Symbolic links are used to bypass the original path where the fileare stored, which contains space characters and causes plink to fail (symbolic links are located in folder simply to make things more tidy). Then we create a file listing all samples id to extract from our newly merged plink files
```bash
run_plink_merge="mkdir temp; 
  cd temp;
  ln -s /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* ./;
  cd ..;
  ls temp/*.bed | sed 's/.bed//g' > plink_files_to_merge.txt;
  awk '{print \$1,\$1}' /mnt/project${dir}random_cohort_participant.tsv | tail -n +2 > iid_list_plink.tsv
  plink --merge-list plink_files_to_merge.txt \
    --keep iid_list_plink.tsv \
    --make-bed --autosome --out ukb22418_c1_22_v2_merged_subset;
  rm plink_files_to_merge.txt;"
```
#### Run command in dxpy:
We use the [Swiss army knife app]() to run our command, since it includes Plink and Regenie softwares. We specify the cloud instance we require  
```bash
dx run swiss-army-knife \
  -icmd="${run_plink_merge}" \
  --tag="plink_merge" \
  --instance-type "mem1_ssd1_v2_x16" \
  --destination="${dir}" --brief --yes
```

############################################################################################################################################################


#### 3 - Plink QC
## Perform standard QC of Plink genotype files

#dir="/Regenie_test/data/"
run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged_subset \
	--autosome --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
	--write-snplist --write-samples --no-id-header --out geno_array_snps_qc_pass"

dx run swiss-army-knife \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.bed" \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.bim" \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.fam" \
	-iin="${dir}random_cohort_participant.tsv" \
	-icmd="${run_plink_qc}" \
	--tag="plink_qc" \
	--instance-type "mem1_ssd1_v2_x16" \
	--destination="${dir}" --brief --yes


#### 4 - Regenie step 1 ####

########################################################################################################################
##### \t not working!!!!!!!!!!!!
##### Is printed as "\t" character and not as tab. Ask Edoardo to check this since you suck with bash

##### Reformat input file to regenie wanted format (can be done maybe at another step?? Directly forst one)
regenie_input_format="echo -e "FID\tIID\tp21022\tp21001_i0" > header.tsv; 
	awk '{print \$1,\$1,\$2,\$3}' random_cohort_participant.tsv | tail -n +2 > temp.tsv;
	cat header.tsv temp.tsv > regenie_pheno_form_input.tsv;
	rm header.tsv temp.tsv"

dx run swiss-army-knife \
	-iin="${dir}random_cohort_participant.tsv" \
	-icmd="${regenie_input_format}" \
	--tag="regenie_input_format" \
	--instance-type "mem1_ssd1_v2_x2" \
	--destination="${dir}" --brief --yes
########################################################################################################################



regenie_step1="regenie --step 1 \
	--bed ukb22418_c1_22_v2_merged_subset \
	--phenoFile regenie_pheno_form_input.tsv \
	--covarFile regenie_pheno_form_input.tsv \
	--extract geno_array_snps_qc_pass.snplist \
	--phenoCol p21001_i0 \
	--covarCol p21022 \
	--out bmi_gwas_test \
	--bsize 1000 --lowmem --qt --loocv --threads 16 --gz"

dx run swiss-army-knife \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.bed" \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.bim" \
	-iin="${dir}ukb22418_c1_22_v2_merged_subset.fam" \
	-iin="${dir}geno_array_snps_qc_pass.snplist" \
	-iin="${dir}regenie_pheno_form_input.tsv" \
	-icmd="${regenie_step1}" \
	--tag="regenie_step1" \
	--instance-type "mem1_ssd1_v2_x16" \
	--destination="/Regenie_test/step1/" --brief --yes

