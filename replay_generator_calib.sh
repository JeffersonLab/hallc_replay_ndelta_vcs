#!/bin/bash
# replay_generator.sh
# Usage: ./replay_generator.sh <run_number>

runnum=$1
nev=2000000
seg=0
jobdir="SLURM/jobs"
logdir="/volatile/hallc/alphaE/ndelta_vcs2/calib/SLURM/logs"

if [ -z "$runnum" ]; then
    echo "Usage: $0 <run_number>"
    exit 1
fi

if (( runnum > 26561 )); then
    prefix="vcs2_production"
else
    prefix="ndelta_production"
fi

mkdir -p ROOTfiles "$jobdir" "$logdir"

rawfile="raw/${prefix}_${runnum}.dat.${seg}"
if [ ! -f "$rawfile" ]; then
    rawfile="cache/${prefix}_${runnum}.dat.${seg}"
fi

if [ ! -f "$rawfile" ]; then
    echo "No raw file found for run $runnum segment $seg"
    exit 1
fi

jobfile="${jobdir}/replay_calib_${runnum}.slurm"

cat > "$jobfile" <<EOF
#!/bin/bash
#SBATCH --account=hallc
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2500
#SBATCH --job-name=replay_${runnum}
#SBATCH --time=24:00:00
#SBATCH --gres=disk:1G
#SBATCH --output=${logdir}/replay_${runnum}.out
#SBATCH --error=${logdir}/replay_${runnum}.err

echo "Running run ${runnum} segment ${seg} on \$(hostname)"

rawfile="raw/${prefix}_${runnum}.dat.${seg}"
if [ ! -f "\$rawfile" ]; then
    rawfile="cache/${prefix}_${runnum}.dat.${seg}"
fi

if [ ! -f "\$rawfile" ]; then
    echo "No raw file found"
    exit 1
fi

cmd="SCRIPTS/COIN/PRODUCTION/replay_production_coin_pElec_hProt.C(${runnum},${nev},${seg})"
hcana -l -q "\$cmd"

outfile="ROOTfiles/coin_replay_production_${runnum}_${nev}_${seg}.root"
if [ ! -f "\$outfile" ]; then
    echo "ERROR: expected output \$outfile not found"
    exit 1
fi
EOF

echo "Created: $jobfile"
echo ""
echo "To submit:"
echo "  sbatch $jobfile"
