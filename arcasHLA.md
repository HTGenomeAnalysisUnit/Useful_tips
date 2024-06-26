# arcasHLA

This is a quick summary of how to run [arcasHLA workflow](https://github.com/RabadanLab/arcasHLA) using a singularity container in our HPC.

## Container image

The singularity image for arcasHLA v0.6.0-0 is available in `/ssu/gassu/singularity/arcas-hla_0.6.0--hdfd78af_0.sif`

## Reference dataset

Reference data folder in the container is `/usr/local/share/arcas-hla-0.6.0-0/dat/`.

The reference data is provided by GAU pre-processed in /ssu/gassu/supporting_files/arcasHLA. You should bind this folder to `/usr/local/share/arcas-hla-0.6.0-0/dat/` inside the container to be able to access the reference data

## Usage 

On a compute node you can run one of the arcasHLA command using the singulatiry container. Remember to bind relevant path in the container, usually

- your woking folder `$PWD`
- the compute node temp folder `/localscratch`
- the arcasHLA reference folder `/ssu/gassu/supporting_data/arcasHLA`

An example of extract command is as follows. Note how the various location are mounted in the container using `-B` argument. The number of threads is set by `-t` and should match the number of CPUs requested for the job

```bash
singularity exec \
    -B /ssu/gassu/supporting_data/arcasHLA:/usr/local/share/arcas-hla-0.6.0-0/dat/ \
    -B $PWD \
    -B /localscratch \
    -B /ssu/gassu \
    /ssu/gassu/singularity/arcas-hla_0.6.0--hdfd78af_0.sif arcasHLA extract \
        --temp $TMPDIR \
        -t 8 -v --log $PWD/extract.log \
        -o $PWD/extract \
        --single --unmapped \
        $PWD/sample.bam
```

## Consideration for single-cell data

Single-cell data are usually aligned as single-end so you should use `--single` option when available. It is also suggested to include `--unmapped` option at the extract stage.

Finally, we may need to aggiust settings for read-length (10X is usually 90bp) and min read support since you may have reduced coverage.

## genotype step example

When running genotype on single-cell data from 10X, the input data is single-end reads with length 90bp. 

It is suggested to first estimate read length and read length SD from the input extracted reads. Here, we estimate these values using only 100k reads to speed-up the process. This is usually fine, but you can remove `head -n 100000` to process the full dataset and get a more accurate estimate. 

You can use a command like this:

```bash
zcat extract/sample.extracted.fq.gz | paste - - - - | head -n 100000 \
    | awk '{sum += length($2); sumsq += (length($2))^2}; END {print "Average:", sum/NR, "; std:", sqrt((sumsq-((sum^2)/NR))/NR)}'

Average: 90 , std: 0
```

If `std` value is zero, you should use a value of one.

Then you can proceed with genotype like this. Note how the read length and standard deviation are provided using `-l` and `-s` options. The number of threads is set by `-t` and should match the number of CPUs requested for the job. Finally, you may need to adjust the minimum N reads support `--min_count`. The default is 70, which is often too high for single-cell data, lowering this number increase sensitivity but also noise..

```bash
singularity exec \
    -B /ssu/gassu/supporting_data/arcasHLA:/usr/local/share/arcas-hla-0.6.0-0/dat/ \
    -B $PWD \
    -B /localscratch \
    /ssu/gassu/singularity/arcas-hla_0.6.0--hdfd78af_0.sif arcasHLA genotype \
        --temp $TMPDIR \
        -t 8 -v \
        --log $PWD/genotype.log \
        -o $PWD/genotype \
        -g A,B,C,DQA1,DQB1,DRB1 \
        --min_count 20 --single -l 90 -s 1 \
        $PWD/extract/sample.extracted.fq.gz
```