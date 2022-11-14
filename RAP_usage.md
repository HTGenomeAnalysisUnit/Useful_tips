# RAP usage

Tips to work in the Research Analysis Platform of UKBB

- [RAP usage](#rap-usage)
  - [Command line access using dxpy](#command-line-access-using-dxpy)
  - [Login and connect to a project](#login-and-connect-to-a-project)
  - [Basics](#basics)
  - [Cloud Workstation App](#cloud-workstation-app)
    - [Start a session](#start-a-session)
    - [Interact with files in the project](#interact-with-files-in-the-project)
    - [Using tools](#using-tools)
    - [Terminate](#terminate)
    - [Logout](#logout)
    - [For further info](#for-further-info)
  - [Run tools with dx](#run-tools-with-dx)
  - [Regenie example](#regenie-example)

## Command line access using dxpy

dxpy is the CLI tool of DNA Nexus that can be used to interact with the platform, like moving files in and out, install tools, run apps etc.

It can be installed easily using conda `conda create -n dxpy_0.327 -c bioconda -c conda-forge dxpy`

For popgen users the tool is available in the conda env `dxpy_0.327` (to activate this conda environment type `conda activate dxpy_0.327`)

## Login and connect to a project

To use dx CLI first login with your RAP credentials using `dx login`.

Then you can see a list of public projects using `dx select --public` or a list of projects you have access to to using `dx select`.

You can connect to a project using `dx select "project_name"`

## Basics

- List files in the project `dx ls`
- Move in a folder `dx cd`
- Upload file into the project `dx upload`
- Configure ssh access for interactive session `dx ssh_login`. Here you can either create a new ssh key or use one that you already have on the system (the tool will prompt you for selection). This is needed before you can use Cloud Workstation App

## Cloud Workstation App

The cloud workstation app allows to start a job on DNANexus and then ssh into the worker so that you obtain an interactive session. You can use it a regular Ubuntu machine, you can install programs, run commands etc. Docker is also available by default.

A useful starting point is in the [official documentation](https://documentation.dnanexus.com/developer/cloud-workstations/cloud-workstation)

Remember to set up you ssh access using `dx ssh_config` before trying to use the CWA.

### Start a session

After you connect to a project, you can start the CWA session using `dx run app-cloud_workstation --ssh`. The tool will prompt you asking which options you want to set. It is better to always set a time limit, while other options can be left unset.

By default, the Cloud Workstation App will launch on a mem1_ssd1_v2_x8 [instance type](https://documentation.dnanexus.com/developer/api/running-analyses/instance-types) which has 8 cores, 16 GB memory, and 180 GB storage. To run the app on a different instance type, use the --instance-type flag for dx run.

`dx run --instance-type mem1_ssd1_v2_x36 --ssh app-cloud_workstation`

Once you get in you need to configure the location to your project using

```bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:
```

### Interact with files in the project

You can copy a file from the parent project to the CWA instance using

`dx download file_name`

And upload a file back into the project using

`dx upload --path "$DX_PROJECT_CONTEXT_ID:"`

Actually, you can download from any project you have access to using the project id like

`dx download project_id:/path/to/file`

### Using tools

In the CWA you can download, compile and run any tool you want. The simplest way to use a tool in the CWA is to have it in a docker container in DockerHub and the you can use

`docker run -v $PWD:/data docker_image you_command /data`

Note that you need to bind the working directory ($PWD) to some path inside the container to be able to see your files.

### Terminate

To terminate a session use `dx terminate $DX_JOB_ID` from inside the session

### Logout

To log out on the command line use `dx logout`

### For further info

See full list of dx commands [here](https://documentation.dnanexus.com/user/helpstrings-of-sdk-command-line-utilities)

## Run tools with dx

In general once a tool is installed it can be run using

```bash
dx run [tool_name] \
  -iifile="filename" #files to be imported from project
  -icmd="[command_to_run]" #String with the command
  --tag="Step1_merge" #A tag for monitorin purposes
  --instance-type "mem1_ssd1_v2_x16" #Instance required for the job
  --destination="/Regenie_test/" #Destination folder within the project to write to
  --brief --yes
```

Keep in mind that for both `iifile` and `destination` you can use absolute path related to the project you are currently logged in or add a project id like `--destination-"project-id:/Data`

In general, files required using `iifile` are visible directly in the working dir of the job. All files from the project are visible in `/mnt/project` read-only.

See [dx run guide](https://documentation.dnanexus.com/user/helpstrings-of-sdk-command-line-utilities#run) for full options

## Regenie example

The following mainly reproduce [this tutorial](https://www.youtube.com/watch?v=762PVlyZJ-U)

1. activate dxpy env `conda activate dxpy_0.327`
2. login `dx login`
3. set up the command

  ```bash
  run_merge="mkdir files;cd files;ln -s /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c* ./;cd ..;ls files/*.bed | sed 's/.bed//g' > files_to_merge.txt;plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out ukb22418_c1_22_v2_merged;rm files_to_merge.txt;"
  ```

  Here we use a trick to avoid issues with spaces that are present in the files path. We first create symlinks in the `files` folder in the working directory and then save the list of files from `files/*.bed` into `files_to_merge.txt`. This file will be used as input for `plink` merge.

4. Run command using dx

  ```bash
  dx run swiss-army-knife \
    -icmd="${run_merge}" \
    --tag="Step1_merge" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="/Regenie_test/" \
    --brief --yes
  ```
  