nohup snakemake --executor slurm --jobs 48 --use-conda --retries 2 --printshellcmds --slurm-no-account > snakemake.log 2>&1 &
