# SRA tips

## Download project data

1. First use pysradb to download the full project metadata including experiment IDs and data URLs

Example:

```bash
pysradb metadata --detailed --expand --saveto PRJNA908079_metadata.csv PRJNA908079
```

2. Using the table produced by the above command the URL to download raw data should be in `ena_fastq_http_1` and `ena_fastq_http_2`. You can use curl to download cycling these columns

Alternatively, you can use ffq with a series of SRR IDs to get the URLs for a specific data type. 

For example, to download the SRA file for a given SRR ID you can use `ffq` and `jq`:

```bash
module load jq
conda activate ffq

ffq --aws SRR22975270 | jq '.[] | select(.filetype == "sra") | .url' | xargs -n 1 curl  
```

## Convert SRA to fastq

To convert SRA to fastq you can use `fasterq-dump` from the SRA toolkit

```bash
module load sratoolkit/3.0.2

fasterq-dump --split-files SRR1234567.sra
```

Note that fasterq-dump do not support direct compression to gz, so you need to compress the output files after downloading
