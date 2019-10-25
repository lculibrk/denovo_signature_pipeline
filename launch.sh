function usage()
{
    echo "launch.sh"
    echo "Launches snakemake workflows on the numbers cluster"
    echo ""
    echo "Usage: bash launch.sh N_CORES SNAKEFILE SMK_ARGS"
    echo ""
    echo -e "N_CORES\tNumber of jobs to submit at a time"
    echo -e "SNAKEFILE\tSnakefile to run"
    echo -e "SMK_ARGS\tOther Snakemake arguments"
}

if [ "$1" == "" ]
then
   usage
   exit
fi

if [ "$HOSTNAME" == "n104" ] || [ "$HOSTNAME" == "n105" ]; then
    echo 'Running on numbers cluster'
    snakemake -p \
        --cluster-config "cluster.json" \
        --drmaa ' --mem-per-cpu={cluster.mem_mb} {cluster.flags} ' \
        --jobs $1 \
        --latency-wait 60 \
        --max-jobs-per-second 5 \
        --restart-times 5 \
        --rerun-incomplete \
        --keep-going \
	--nolock \
	-s $2 \
	${@:3}
    mkdir -p logs/cluster_logs
    mv slurm*.out logs/cluster_logs

elif [ "$HOSTNAME" == "login-apollo.hpc.bcgsc.ca" ]; then
    echo "Running on APOLLO cluster"
    mkdir -p logs/sge_logs_error logs/sge_logs_output
    snakemake -p \
        --cluster-config "config/cluster.json" \
        --cluster "qsub -q arc.q -P arc.prj -V -N msig_timing -e logs/sge_logs_error -o logs/sge_logs_output -l h_vmem={cluster.mem_mb}M" \
        --jobs 500 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete

else
    echo "Running in non-cluster mode"

    snakemake -p \
        --cores $1 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete \
        -s $2
fi
