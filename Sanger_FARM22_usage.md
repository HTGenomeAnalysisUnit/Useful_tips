# Connect and use FARM22 cluster

To be able to use the Sanger cluster you need
- Sanger user name
- Sanger user password
- Be added to the cardinal_analysis group in bsub (this may need a request to the service desk and you have to take the FARM training course)

## Table of contents

- [Connect and use FARM22 cluster](#connect-and-use-farm22-cluster)
  - [Table of contents](#table-of-contents)
  - [Connect to FARM22 cluster](#connect-to-farm22-cluster)
    - [Use Teleport on the web](#use-teleport-on-the-web)
    - [Set up teleport on your computer to connect using SSH](#set-up-teleport-on-your-computer-to-connect-using-ssh)
    - [Access Sanger server in SSH](#access-sanger-server-in-ssh)
    - [Set up a proxy to access internal Sanger web pages](#set-up-a-proxy-to-access-internal-sanger-web-pages)
  - [Set up connection in Windows](#set-up-connection-in-windows)
  - [Copy data from FARM22 cluster](#copy-data-from-farm22-cluster)
  - [Use VS Code on FARM22 cluster](#use-vs-code-on-farm22-cluster)
  - [Configure your environment on FARM22](#configure-your-environment-on-farm22)
    - [Access basic software](#access-basic-software)
    - [Add new software](#add-new-software)
    - [Use jupyter notebook](#use-jupyter-notebook)
    - [Basics of the scheduler](#basics-of-the-scheduler)
      - [Main scheduler commands](#main-scheduler-commands)
      - [interactive jon template](#interactive-jon-template)
      - [batch job script template](#batch-job-script-template)
  - [CARDINAL data](#cardinal-data)
    - [Summary of data organization](#summary-of-data-organization)
    - [CARDINAL dashboard](#cardinal-dashboard)

## Connect to FARM22 cluster

Sanger's access system is now based on Teleport. There are multiple ways to access using either a browser interface or the command line.

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) for details on how to set-up teleport and ssh connection for the new system.

### Use Teleport on the web

The first way is to log in with your Sanger credential at `https://sanger.okta.com/` and then select the Teleport application. After login again with your user credentials and Okta authentication, you can select to open a terminal on either the gen22 or farm22 system. 

You can then work freely using the command line directly on your browser.

See this [video tutorial](https://app.tango.us/app/workflow/Connect-to-the-Sanger-server-using-Teleport-on-the-web-ea3b4c970ee4445f915bc2226ffa400c) for a demonstration of how to use Teleport on the web.

### Set up teleport on your computer to connect using SSH

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) for details on how to set up teleport and ssh connection from your laptop.

Here is a quick start summary of the main steps.

1. Install Teleport on your computer using the appropriate installer from [Teleport download page]( https://goteleport.com/download/#install-links)
2. Authenticate your user with Okta with this command

   `tsh login --proxy=portal.sanger.ac.uk:443 --auth=okta`
   
4. Generate a config file with this command `tsh config >> ~/.ssh/config`
5. Note down the location of the Identity and Certificate file that you can get using
    - Identity file: `grep IdentityFile ~/.ssh/config | grep portal | head -n1`
    - Certificate file: `grep CertificateFile ~/.ssh/config | grep portal | head -n1`
6. Create a folder to store control paths using `mkdir ~/.ssh/controlpaths`
7. Open the file `~/.ssh/config` in a text editor and add the following at the beginning

   ```
   ControlMaster auto
   ControlPath ~/.ssh/controlpaths/s_%C
   TCPKeepAlive yes
   AddKeysToAgent yes
   ForwardAgent yes
   ServerAliveInterval 20
   ServerAliveCountMax 5
   ConnectTimeout 45
   ```  

9. In the file `~/.ssh/config` add the following block after the line `# End generated Teleport configuration`

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

Host farm22
    User <YOUR SANGER USER NAME HERE>
    HostName farm.internal.sanger.ac.uk
    ProxyJump jammy-dev64
```

### Access Sanger server in SSH

Once you have completed the configuration you can log in with Okta using

`tsh login --proxy=portal.sanger.ac.uk:443 --auth=okta`

Then you can connect to Sanger FARM22 from the command line using

`ssh farm22`

**NB.** The okta token expires every 8h so you must log in again after this time interval.

### Set up a proxy to access internal Sanger web pages

Once you have completed the configuration above 

1. Open a tunnel to the Sanger server using the configured SSH connection. You can use a command like this that will open a connection in the background on your computer

   `ssh -Nf jammy-dev64`

   You can close the tunnel using this command

   `ssh -O exit jammy-dev64`
   
3. In Firefox, go to Settings, search for `proxy` and then select Network settings and set your proxy as in this image

![Proxy settings](<Sanger_FARM22_proxy_settings.png>)

Now you should be able to browse the Sanger internal web pages.

**NB.** It is suggested to change the proxy in a dedicated browser different from your default one. Once the proxy is set, all traffic will be redirected through the proxy so the browser will may not be able to connect to normal Internet websites.

## Set up connection in Windows

After configuring Teleport access as described above, you can access to the Sanger FARM22 using this teleport command:

`tsh login --proxy=portal.sanger.ac.uk:443 --auth=okta`

To configure the tunnelling in Windows and being able to access Sanger internla resources, you need to configure a tunnel in [Putty](https://www.putty.org/). You can follow these steps:

1. Download putty [last version 1.8](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
2. `tsh puttyconfig YOUR_USERNAME@farm22-head2`
3. Open putty
4. In the Putty configuration, select the newly configured proxy for your session and click load

![screenshot](<putty_one.png>)

5. Go to Connection > SSH > Tunnels: insert `L3128` in Source port and `webcache.sanger.ac.uk:3128` in Destination and click add

![screenshot](<putty_two.png>)

6. go back to Session, insert farm22-head2 into the Host Name and save.

![screenshot](<putty_three.png>)

Now next time you need to access to Sanger HPC and open the tunnel follow these steps:

1. open putty
2. click on farm22-head2 (proxy:portal.sanger.ac.uk) in the saved session 
3. load
4. open
5. you should see a browser opening and Teleport will ask for Okta authentication.
6. you will be able to Navigate to Sanger webpages in Mozilla if you configured as explained above. 

## Copy data from FARM22 cluster

Once you completed the configuration above you can use Teleport to copy files from FARM22 through scp. 

Example of transferring a file from FARM22 to your local computer

```bash
tsh scp <username>@farm22-head1:/my/path/test.txt ./
```

Or the other way round from your local system to FARM22

```bash
tsh scp ./test.txt <username>@farm22-head1:/my/path
```


## Use VS Code on FARM22 cluster

If you performed the configuration described above, when using the VSCode SSH remote extension, you should see a `farm22` host that you can connect to. Using this you will open a remote VSCode session on farm22 login node.

**NB.** Since you are connected to the login node, you can use this session to edit your code, but you MUST NOT perform any heavy computation in this session (including running intensive analysis in a jupyter notebook).

See the [Sanger documentation](https://sanger.freshservice.com/support/solutions/articles/53000059623) on how to configure VSCode to connect to FARM servers.

## Configure your environment on FARM22

### Access basic software

Many tools are provided directly by HGI and other groups as modules. Check the list of available tools using `module avail` or `module avail <tool name>`.

To enable the use of conda you need to load the relevant module

```bash
module load HGI/common/conda/module
```

Then you can use the standard conda commands. 

**.NB** To activate a conda environment one must use the `source` command like in this example 

```bash
source activate /full/path/to/conda/env
```

### Add new software

- Additional software must be installed in `/software/cardinal_analysis/ht`.
- In this location, we have created a `conda_envs` sub-folder that can be used when creating new conda environments
- Similarly, we also created a `singularity` subfolder to store singularity containers

### Use jupyter notebook

Jupyter notebooks should preferentially be run using the Sanger Jupyter Hub. Once you set up the proxy connection as described above, you can access the Jupyter Hub at [https://jupyter.internal.sanger.ac.uk/](https://jupyter.internal.sanger.ac.uk/).


### Basics of the scheduler

Sanger uses bsub scheduler.

#### Main scheduler commands

Here we listed the main scheduler commands with the main options

**bsub** (used to submit jobs)

- `-e / -o`: STDERR and STDOUT log files (you can use %J to add job id to the file name and %I to add the index of the array job)
- `-Is`: interactive job
- `-J`: job name
- `-n`: n cores
- `-M`: memory
- `-q`: queue
- `-R`: resources 
- `-G`: group (in theory we all part of cardinal_analysis group)

**bjobs** (list jobs)

- `-p`: reason for pending
- `-l`: long details on the job

**bqueues** (show queues status, use this to see available queues)

- `-w` status info
- `-l` jobid: full stats

**bhosts: list hosts**

- `-w` more info

**bugroup** (show LSF bsub groups)

To see which LSF groups you are member of you can use: `bugroup | grep $USER`

#### interactive job template

To request an interactive session you use a command like this. Adjust the memory and cpu requirements as needed.

```
bsub -Is -n 1 -M 4G -q normal -G cardinal_analysis -R "select[mem>4G] rusage[mem=4G] span[hosts=1]" -J test /bin/bash
```

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
- `analysis`: contains the analysis results for each user. **Run all your analyses here.**
- `qc`: main data organized by tranche
- `freezes`: static copy of the data freezes as discussed in meetings. Most of the time you want to use one of these releases.

### Summary of data organization

1. Each tranche has several donor pools that are stored in the `Donor_Quantification` subfolder in a tranche folder. Like: 

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

There are lots of plots displaying how the deconvolution, cell type assignment etc. in the [official CARDINAL dashboard](https://apps.hgi.sanger.ac.uk/scrna/).

Choose the ‘cardinal analysis’ tab to see the QC plots for each stage of the analysis for each tranche:
