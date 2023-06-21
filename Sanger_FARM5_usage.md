# Connect and use FARM5 cluster

To be able to user Sanger cluster you need
- Sanger user name
- Sanger SSH password
- Sanger user password (note this is different from the above)
- Be added to the cardinal_analysis group in bsub (this may need a request to service desk and you have to take the FARM5 training course)

## Connect to FARM5 cluster

The address of sanger cluster access point from external network is `ssh.sanger.ac.uk`.

You can connect with usual ssh command, using your SSH password:

```bash
ssh username@ssh.sanger.ac.uk
```

Then yo will be prompted to select the server location. Here you can use `farm5` and then use you Sanger user password to connect to the cluster.

### Set up proxy to access internal sanger web pages

If you need to access internal Snager web pages (like Matiss's app or some confluence documentation) you need to set up a tunnel and then use a proxy on your browser (I suggest to use Firefox or anyway a different browser from your usual one so you don't lose access to normal Internet service on your default browser).

First, connect to the cluster, but using following command:

```bash
ssh -L 3128:webcache.sanger.ac.uk:3128 username@ssh.sanger.ac.uk
```

Then on Firefox, goes to Settings, search for `proxy` and then select Network settings and set your proxy as in this image

![Proxy settings](<Sanger_FARM5_proxy_settings.png>)

Now you should be able to browse Sanger internal web pages.

## Copy data from FARM5 cluster

Copying of files needs to be done through SSH port forwarding (a tunnel). 

This is setup as follows, from the client machine:

```bash
ssh -L 2222:internal_sanger_host:22 username@ssh.sanger.ac.uk 
```

Under normal circumstances, the `internal_sanger_host` will be the same host you specify when prompted after your initial login, thus usally `farm5-head1`.

This then allows you to copy files to the Sanger machine (or from the Sanger machine) using `scp` or `scftp`. When using the copy commands, you need to connect to 'localhost' (usually IP 127.0.0.1) on port '2222'. For example:

```bash
#Copy a file to Sanger
scp -P2222 filename username@localhost:filename

#Copy a file from Sanger
scp -P2222 username@localhost:filename filename
```

## Use VS Code on FARM5 cluster

If you have the SSH connection extension in VS code, you can work remotely in FARM5 following these steps.

First, add this to the `~/.ssh/config` file in your computer

```bash
Host farm5
  HostName 127.0.0.1
  Port 5022
  User username
```

Then create a tunnel to FARM5 cluster as described above, using port 5022 (you can change this to any port you like, but then you need to change the port in the `~/.ssh/config` file as well)

```bash
ssh -L ssh -L 5022:farm5-head2:22 eg20@ssh.sanger.ac.uk
```

Then reload VScode and you should see a `farm5` host in the available ssh remotes and you can connect using your Sanger user password.

## Configure your environment

### Access basic software

This is needed to access basic tools like conda, and some basic modules.

First do this

```bash
echo "hgi" >> ~/.softwarerc
```

Then add this to your `.bashrc`

```bash
module use --append /software/modules
source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
```

See also https://confluence.sanger.ac.uk/pages/viewpage.action?spaceKey=HGI&title=Software+on+the+Farm 

**NB.** Additional softwares must be installed in /software/cardinal_analysis/

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

Main location of CARDINAL data is: `/lustre/scratch123/hgi/projects/cardinal_analysis/` qc

Inside this directory there are several subdirectories. Main ones are:
- `analysis`: contains the results of the analysis for each user. **Run all your analysed here.**
- `qc`: main data organized by tranche
- `freezes`: static copy of the data freezes as dicussed in meetings

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