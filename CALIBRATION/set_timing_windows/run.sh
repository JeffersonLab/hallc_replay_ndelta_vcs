#This PR contains the hodoscope timing toolkits.
#
#After this PR is updated, 
#
## sync the forked repository
#
#```
#git fetch upstream
#git checkout main
#git merge upstream/main
#```
#
## check calib_run_list_by_run_period.csv file for which run is assigned to each student.
#
#- Penn: 26369
#- Chathuri: 26382
#- Thareendra: 26549
#- Nazmus: 26567
#- Suman: 27064
#
## run the timing window script
#```
#cd CALIBRATION/set_timing_windows/
#root -l 'hodo_timediff_cut_app.C([run_number], 26772)'
#```
#
#to produce the json files. 
#After the task is done, execute the following
#```
#root -l -b -q 'hodo_timediff_cut_app.C([run_number], [run_number], true)'
#```
#
## send a PR
#
#Include the json file and **pdf** in the PR!!
#
#You can use
#
#```
#git add -f [run_number]_all_summary.pdf
#```

# After the time window was surveyed, 
# ```
# python3 merge_hodo_qa.py
# ```
# was used to generate the related hodoscope cuts.
# Other related parts were coped from the existing hhodo_cuts.param, phodo_cuts.param



# Similarly, root -l -b -q "other_det_timediff_cut_app.C(\"all\", ${runnum})";
# will open the other detector timing cut window
# However, this script was not used.
# Existing time window parameters were good enough.

# Below is the verification attempt in batch mode, but we ended up using slurm to submit jobs.

#for runnum in 26088 26169 26772; do     root -l -b -q "hodo_timediff_cut_app.C($runnum, $runnum, true)";   done
#for runnum in $(ls ../../ROOTfiles/ | grep -oP 'coin_replay_production_\K[0-9]+' | sort -un); 
#  do
#  root -l -b -q "other_det_timediff_cut_app.C(\"all\", ${runnum}, -1, true)";
#  root -l -b -q "hodo_timediff_cut_app.C(${runnum}, 0, true)";
#done
#echo "Do not runnum this... too slow. do verification_generator.sh"

