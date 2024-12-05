# RAP usage

Tips to work in the Research Analysis Platform of UKBB

- [RAP usage](#rap-usage)
  - [Command line access using dxpy](#command-line-access-using-dxpy)
  - [Before you start - Login and connect to a project](#before-you-start---login-and-connect-to-a-project)
  - [Navigate a project files](#navigate-a-project-files)
  - [Run tools with dx](#run-tools-with-dx)
  - [The Cloud Workstation App](#the-cloud-workstation-app)
    - [Start a session](#start-a-session)
    - [Interact with files in the project](#interact-with-files-in-the-project)
    - [Use tools](#use-tools)
    - [Terminate](#terminate)
    - [Further info](#further-info)

## Command line access using dxpy

`dxpy` is the CLI tool of DNA Nexus that can be used to interact with the platform, like moving files in and out, install tools, run apps etc.

It can be installed easily using conda `conda create -n dxpy_env -c bioconda -c conda-forge dxpy`

## Before you start - Login and connect to a project

To use dx CLI first you need to login with your RAP credentials using `dx login`.

You also need to configure ssh keys to access interactive session: `dx ssh_login`. Here you can either create a new ssh key or use one that you already have on the system (the tool will prompt you for selection). This is needed before you can use Cloud Workstation App.

To see the list of available projects for your user and connect to a specific one you can use 

- `dx select --public` : list publicly available projects
- `dx select --level VIEW` : list all projects to you are a member of
- `dx select` : list projects you have contributing access to

You can connect to a project using `dx select "project_name"` or selecting one from the list after one of the above commands

## Navigate a project files

The `dx` tool can be used to navigate project folder, list files and basic file operations like move, copy, delete. The syntax use the standard linux commands like `ls`, `cd`, `mv`, `cp`, `rm` etc.

Keep in mind that your space depends on the project you are connected to.

Examples:

- List files in the project `dx ls`
- Move in a folder `dx cd /folder/in/the/project`, also `dx cd ..` works
- Upload file into the project `dx upload`
- See the current directory of the project `dx pwd`

To search files in the project you can use `dx find data`. This works with some specific options, the most relevant to know are `--name` to set a pattern or name to search and `--class` to set the type of object to search for (i.e. applet, file, workflow, ...). Example command to search for files: `dx find data --name "*.bim" --class file`

## Run tools with dx

In general once a tool is installed it can be run using

```bash
dx run [tool_name] \
  -iin="filename" #files to be imported from project
  -icmd="[command_to_run]" #String with the command
  --tag="Step1_merge" #A tag for monitorin purposes
  --instance-type "mem1_ssd1_v2_x16" #Instance required for the job
  --priority "normal" #Priority of the job, you can use high, normal, low. High cost more
  --destination="/Regenie_test/" #Destination folder within the project to write to
  --brief --yes
```

Keep in mind that for both `iin` and `destination` you can use absolute path related to the project you are currently logged in or add a project id like `--destination-"project-id:/Data`

In general, files required using `iin` are visible directly in the working dir of the job. All files from the project are visible in `/mnt/project` read-only.

See [dx run guide](https://documentation.dnanexus.com/user/helpstrings-of-sdk-command-line-utilities#run) for full options

## The Cloud Workstation App

The cloud workstation app allows to start a job on DNANexus and then ssh into the worker so that you obtain an interactive session. You can use it a regular Ubuntu machine, you can install programs, run commands etc. Docker is also available by default.

A useful starting point is in the [official documentation](https://documentation.dnanexus.com/developer/cloud-workstations/cloud-workstation)

Remember to set up you ssh access using `dx ssh_config` before trying to use the CWA.

**Interactive session can quickly become expensive on large instances**. Use them wisely and mainly for quick data exploration or debugging.

**All files generated during a session are lost unless manually uploaded to you project**. Remember to upload any relevant data before terminating the session.

### Start a session

After you connect to a project, you can start the CWA session using `dx run app-cloud_workstation`. The tool will prompt you asking which options you want to set. **It is better to always set a time limit**, this will avoid to incurr in excessive costs in case you forgot to shutdown the session.

By default, the Cloud Workstation App will launch on a mem1_ssd1_v2_x8 [instance type](https://documentation.dnanexus.com/developer/api/running-analyses/instance-types) which has 8 cores, 16 GB memory, and 180 GB storage. To run the app on a different instance type, use the `--instance-type` flag for dx run.

```bash
dx run --instance-type mem1_ssd1_v2_x2 app-cloud_workstation`
```

This command will show a job ID in the form `Job ID: job-Gx8ppgQJ08pf6FYK9Yj6fP8q` . You can then use this id to ssh into the running machine and open your terminal

```bash
dx ssh job-Gx8ppgQJ08pf6FYK9Yj6fP8q
```

Once you get in you need to configure the location to your project using:

```bash
unset DX_WORKSPACE_ID
dx cd $DX_PROJECT_CONTEXT_ID:
```

### Interact with files in the project

In the interactive virtual environment, the files from your project are not immediately accessible. To access your porject files you must use `dx` to navigate to the project folder of interest and then use `dx download` and `dx upload` to move files in and out of the CWA instance.

You can copy a file from the project data to the CWA instance using

```bash
dx download file_name #Copy file from the current working directory in the project
dx download /path/to/file #Copy file from a specific path in the project

#Actually, you can download from any project you have access to using the project id like
dx download project_id:/path/to/file
```

Similarly, you can upload a file back into the project using

```bash
#Copy a file in the current local folder to the current working directory in the project
dx upload local_file_name 

#Copy a file to a specific path in the project. Note that you have to provide full path including file name
dx upload --path "$DX_PROJECT_CONTEXT_ID:/my/file" local_file_name 
```

### Use tools

In the CWA you can download, compile and run any tool you want. The simplest way to use a tool in the CWA is to have it in a docker container in DockerHub and the you can use

```bash
docker run -v $PWD docker_image your_command /data
```

Note that you need to bind the working directory ($PWD) to some path inside the container to be able to see your files.

### Terminate

To terminate a session use `dx terminate $DX_JOB_ID` from inside the session

### Further info

See full list of dx commands [here](https://documentation.dnanexus.com/user/helpstrings-of-sdk-command-line-utilities)
  
