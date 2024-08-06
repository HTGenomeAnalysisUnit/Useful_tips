# Connect and use FARM5 cluster

To be able to user Sanger cluster you need
- Sanger user name
- Sanger SSH password
- Sanger user password (note this is different from the above)
- Be added to the cardinal_analysis group in bsub (this may need a request to service desk and you have to take the FARM5 training course)

## Connect to FARM22 cluster

Sanger access system is now based on Teleport. There are multiple ways to access using either browser interface or command-line.

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) for details on how to set-up teleport and ssh connection for the new system.

### Use Teleport on the web

The first way is to login with your Sanger credential at `https://sanger.okta.com/` and then select Teleport application. After log-in afgain with your user credentials and Okta authentication, you can select to open a terminal on either gen22 or farm22 system. 

You can then work freely using the command-line directly on your browser.

See this [video tutorial](https://app.tango.us/app/workflow/Connect-to-the-Sanger-server-using-Teleport-on-the-web-ea3b4c970ee4445f915bc2226ffa400c) for a demonstration of how to use Teleport on the web.

### Set up teleport on you computer to connect using SSH

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) for details on how to set-up teleport and ssh connection from you laptop.

Here is a quick start summary of the main steps.

1. Install Teleport on your computer using one of appropriate installer from [Teleport download page]( https://goteleport.com/download/#install-links)
2. Authenticate your user with Okta with this command: `tsh login --proxy=portal.sanger.ac.uk:443 --auth=okta`
3. Generate a config file with this command `tsh config > ~/.ssh/teleport_config`
4. Note down the location of the Identity and Certificate file that you can get using grep `IdentityFile ~/.ssh/teleport_config | grep portal | head -n1` and `grep CertificateFile ~/.ssh/teleport_config | grep portal | head -n1`, respectively.
5. Open the file `~/.ssh/teleport_config` in a text editor and add the following after the line `# End generated Teleport configuration`

```
Host jammy-dev64 focal-dev64
    HostName %h
    Port 3022
    UserKnownHostsFile ~/.ssh/teleport_known_hosts
    LocalForward 3128 webcache.sanger.ac.uk:3128
    ForwardX11 yes
    ForwardAgent yes
    User <YOUR SANGER USER NAME HERE>
    ProxyCommand "/usr/local/bin/tsh" proxy ssh --cluster=portal.sanger.ac.uk --proxy=portal.sanger.ac.uk:443 %r@%h:%p
    IdentityFile "<YOUR IDENTITY FILE PATH HERE>"
    CertificateFile "<YOUR CERTIFICATE FILE PATH HERE>"
    HostKeyAlgorithms +ssh-rsa*,rsa-sha2-512
    PubkeyAcceptedKeyTypes +ssh-rsa*
    DynamicForward 8888
```

### Access Sanger server in SSH

Once you completed the configuration you can login with Okta using

`tsh login --proxy=portal.sanger.ac.uk:443 --auth=okta`

And then you can connect to Sanger FARM22 from the command line using

`tsh ssh eg20@farm22-head1`

**NB.** The okta token expires every 8h so you will have to login again after this time interval.

### Set up proxy to access internal sanger web pages

Once you completed the configuration above 

1. Connect to Sanger server using the configured SSH connections
2. In Firefox, goes to Settings, search for `proxy` and then select Network settings and set your proxy as in this image

![Proxy settings](<Sanger_FARM22_proxy_settings.png>)

Now you should be able to browse Sanger internal web pages.

**NB.** It is suggested to change the proxy in a dedicated browser different from the one you use usually. Once the proxy is set, all traffic will be redirected through the proxy so the browser will may not be able to connect to normal Internet websites.

## Copy data from FARM22 cluster

Once you completed the configuration above you can use Teleport to copy files from FARM22 thrugh scp. 

Example of transferring a file from FARM22 to local

```bash
tsh scp <username@farm22-head1:/my/path/test.txt ./
```

Or the other way round from local to FARM22

```bash
tsh scp ./test.txt <username@farm22-head1:/my/path
```


## Use VS Code on FARM22 cluster

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) on how to configure VSCode to connect to FARM servers.

## Configure your environment on FARM22

### Access basic software

Many tools are provided directly by HGI and other groups as modules. Check the list of available tools using `module avail` or `module avail <tool name>`.

To enable the use of conda you need to load the relevant module

```bash
module load HGI/common/conda/module
```

Then you can use the standard conda commands. 

**.NB** To activate a conda environment one must use `source` command like in this example 

```bash
source acitvate /full/path/to/conda/env
```

### Add new software

- Additional software must be installed in `/software/cardinal_analysis/ht`.
- In this location we have created a `conda_envs` sub-folder that can be used when creating new conda environments
- Similarly we also created a `singularity` subfolder to store singulatiry containers

### Use jupyter notebook

Jupyter notebook should always be run using the Sanger Jupyter Hub. Once you set up the proxy connection as described above, you can access the Jupyter Hub at [https://jupyter.internal.sanger.ac.uk/](https://jupyter.internal.sanger.ac.uk/).


### Basics of the scheduler

Sanger uses bsub scheduler.

#### Main scheduler commands

Here we listed the main scheduler commands with main options

**bsub** (used to submit jobs)
-e / -o: STDERR and STDOUT log files (you can use %J to add job id to the file name and %I to add the index of the array job)
-I: interactive job
-J: job name
-n: n cores
-M: memory
-q: queue
-R: resources 
-G: group (in theory we all part of cardinal_analysis group)

**bjobs** (list jobs)
-p: reason for pending
-l: long details on the job

**bqueues** (show queues status, use this to see available queues)
-w status info
-l jobid: full stats

**bhosts: list hosts**
-w more info

**show_my_lsf_groups** (show your bsub groups)

#### batch job script template

```bash
#!/bin/bash 
#BSUB -n 4
#BSUB -M 256G
#BSUB -q normal
#BSUB -G cardinal_analysis
#BSUB -R "select[mem>256G] rusage[mem=256G] span[hosts=1]"
#BSUB -o batch_analysis.%J.log
#BSUB -e batch_analysis.%J.log
#BSUB -J batch_analysis

#Add the following line if you need to activate conda envs in your script
source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
```

## CARDINAL data

Main location of CARDINAL data is: `/lustre/scratch123/hgi/projects/cardinal_analysis/`

Inside this directory there are several subdirectories. Main ones are:
- `analysis`: contains the results of the analysis for each user. **Run all your analysed here.**
- `qc`: main data organized by tranche
- `freezes`: static copy of the data freezes as discussed in meetings. Most of the time you want to use one of these releases.

### Summary of data organization

1. Each tranche has a number of donor pools that are stored in the `Donor_Quantification` subfolder in a tranche folder. Like: 

```bash
qc/Cardinal_46112_Nov_02_2022/Donor_Quantification/
├── CRD_CMB13219750
├── CRD_CMB13219751
├── CRD_CMB13219752
├── CRD_CMB13219753
├── CRD_CMB13219754
├── CRD_CMB13219755
├── CRD_CMB13219756
├── CRD_CMB13219757
├── CRD_CMB13219758
├── CRD_CMB13219759
└── CRD_CMB13219760
```

And each donor pool has been deconvoluted into donors. Individual h5ad files are stored in each pool subfolder. 

The link between labels like ‘donor0’ and the actual genotype labels for ELGH etc is in the file at the tranche level. Like `/lustre/scratch123/hgi/projects/cardinal_analysis/qc/ELGH_26thMay_2022/Donor_Quantification/assignments_all_pools.tsv`


### CARDINAL dashboard

There are lots of plots displaying how the deconv, cell type assignment etc in the [official CARDINAL dashboard](https://apps.hgi.sanger.ac.uk/scrna/).

Choose the ‘cardinal analysis’ tab to see the QC plots for each stage of the analysis for each tranche: