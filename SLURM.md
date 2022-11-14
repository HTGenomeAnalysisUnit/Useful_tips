# SLURM

## Monitor server load

1. Get summary of status per node, reporting name, free mem, tot mem, state, N CPUs, CPUs state (Allocated/Idle/Other/Total)

```bash
sinfo -o "%n %e %m %a %c %C"
```

2. Get info on the partitions

```bash
sinfo -l
```

3. List active/queued jobs

```bash
squeue #all users

squeue -u username #for a specific user
```

## Details on a job

This returns full details on a queued/running job given the job ID

```bash
scontrol show jobid -dd [job_number]
```
