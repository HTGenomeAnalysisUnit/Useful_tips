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

## List active/queued jobs

```bash
squeue #all users

squeue -u username #for a specific user
```

## Details on a job

This returns full details on a queued/running job given the job ID

```bash
scontrol show jobid -dd <job_number>
```

## Monitor resources

When a job is finished you can inspect the resource usage and the efficency using `seff`. This works also when the job is still running, but the reported stats may be not reliable

```bash
seff <job_number>
```

If you want to get real-time stats for a running job you can use `sstat` instead

```bash
sstat <job_number>

#To return an easy to parse format with columns separated by |
sstat --parsable <job_number>
```

Eventually, you can set which information to return using `--format` and comma separated list of fields. Available fields can be listed using `sstat --helpformat`.

## Cancel a job based on jobname

For a single job name

```bash
scancel -n <jobname>
```

Or to capture all jobs with a pattern

```bash
squeue -u $USER -o "%.10i %.50j" | grep nf-CALL | awk '{print $1}' | xargs scancel
```
