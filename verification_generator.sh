#!/bin/bash
# verification_generator.sh
# Usage: ./verification_generator.sh <run_number>

# foreach file ( `ls CALIBRATION/set_reftimes/reftime_qa/tcoin_*` )
# ./verification_generator.sh `echo $file | sed 's/CALIBRATION\/set_reftimes\/reftime_qa\/tcoin_//g;s/.param//g'`
# sbatch SLURM/jobs/verification_`echo $file | sed 's/CALIBRATION\/set_reftimes\/reftime_qa\/tcoin_//g;s/.param//g'`.slurm
# end

runnum=$1
jobdir="SLURM/jobs"
logdir="/volatile/hallc/alphaE/ndelta_vcs2/calib/SLURM/logs"
jobfile="${jobdir}/verification_${runnum}.slurm"



cat > "$jobfile" <<EOF
#!/bin/bash
#SBATCH --account=hallc
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2500
#SBATCH --job-name=verification_${runnum}
#SBATCH --time=24:00:00
#SBATCH --gres=disk:1G
#SBATCH --output=${logdir}/verification_calib_${runnum}.out
#SBATCH --error=${logdir}/verification_calib_${runnum}.err

cd CALIBRATION/set_reftimes
echo "Running run ${runnum} on \$(hostname)"

root -l -b -q 'reftime_cut_app.C("all",${runnum},${runnum},true)'
EOF

echo "Created: $jobfile"
echo ""
echo "To submit:"
echo "sbatch $jobfile"
