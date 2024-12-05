# Use applet from command line

- [Use applet from command line](#use-applet-from-command-line)
  - [General principles](#general-principles)
  - [Access to files](#access-to-files)
  - [The swiss army knife](#the-swiss-army-knife)
    - [Regenie example](#regenie-example)

## General principles

You can run any of the tools (applets) installed in your project from the command line using the `dx run` command. The general syntax is:

```bash
dx run [applet_name] \
  -iin="filename" #files to be imported from project into the applet before run (aka input files)
  --tag="tag_name" #tag to identify the run
  --instance-type "mem1_ssd1_v2_x16" #instance type to use
  --priority "normal" #Priority of the job, you can use high, normal, low. High cost more
  --destination="/path/to/folder/" #full path to a folder in the project to store generated data
  --brief --yes #flags to avoid confirmation prompts
```

The `-i` option accepts any pair of `key=value` to configure the inputs needed by the specific applet, as defined in the applet documentation you can read in the RAP interface.

Alternatively, you can define all the inputs in a JSON file and use the `-f JSON_FILE` to configure the inputs.

## Access to files

To access files when running the applet you can specify the needed files using `-iin`, this is usually the preferred way. 

Otherwise, data will be available read-only from `/mnt/project` in the applet container. 

All new files created during the processing will be copied to the destination folder in the project after successful completion of the applet.

## The swiss army knife

Swiss army knife is the main tool available in the RAP, and it contains many useful tools (bcftools, samtools, regenie, bgzip, tabix, etc.). You can use it to run any of these tools from the command line. It can also be used to run custom docker images.

```bash
dx run swiss-army-knife \
  -iin="filename" #files to be imported from project into the applet before run (aka input files). Can be repeated
  -icmd="[command_to_run]" #String with the command
  --tag="tag_name" #tag to identify the run
  --instance-type "mem1_ssd1_v2_x16" #instance type to use
  --destination="/path/to/folder/" #full path to a folder in the project to store generated data
  --brief --yes #flags to avoid confirmation prompts
```

You can add `-iimage docker_image` to run a custom docker image. The image should be available in a public repository or in the DNA Nexus platform. In this case commands provided by `-icmd` will be run inside the docker container.

```bash
dx run swiss-army-knife \
  -iin="filename" #files to be imported from project into the applet before run (aka input files). Can be repeated
  -icmd="[command_to_run]" #String with the command]
  -iimage "docker_image" #docker image to use
  --tag="tag_name" #tag to identify the run
  --instance-type "mem1_ssd1_v2_x16" #instance type to use
  --priority "normal" #Priority of the job, you can use high, normal, low. High cost more
  --destination="/path/to/folder/" #full path to a folder in the project to store generated data
  --brief --yes #flags to avoid confirmation prompts
```

### Regenie example

The following mainly reproduce [this tutorial](https://www.youtube.com/watch?v=762PVlyZJ-U)

1. activate dxpy env `conda activate dxpy_0.327`
2. set up the command

  ```bash
  run_merge="mkdir files;cd files;ln -s /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c* ./;cd ..;ls files/*.bed | sed 's/.bed//g' > files_to_merge.txt;plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out ukb22418_c1_22_v2_merged;rm files_to_merge.txt;"
  ```

  Here we use a trick to avoid issues with spaces that are present in the files path. We first create symlinks in the `files` folder in the working directory and then save the list of files from `files/*.bed` into `files_to_merge.txt`. This file will be used as input for `plink` merge.

3. Run command using dx

  ```bash
  dx run swiss-army-knife \
    -icmd="${run_merge}" \
    --tag="Step1_merge" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="/Regenie_test/" \
    --brief --yes
  ```