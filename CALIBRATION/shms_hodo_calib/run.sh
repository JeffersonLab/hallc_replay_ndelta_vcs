# Each app supports interactive/ noninteractive mode, but only noninteractive mode was used for TW and lambda part

# TW calibration
# root -l -b -q 'timewalk_calib_app.C(0, 0, "vanilla", true, "26483,26484,26485,26486,26487,26488", 10.0, "", "./timewalk_qa", "../../PARAM", 1200.0, 20.0, 300.0, true)'
# compareOnly is now true, so it will not update the TW fit. compareOnly=kFALSE can be used as an argument, or one can update timewalk_calib_app.C
# Or, to update the TW calibration one can do simply,
# root -l -b -q 'timewalk_calib_app.C(0, 0, "vanilla", true, "26483,26484,26485,26486,26487,26488", 10.0, "", "./timewalk_qa", "../../PARAM", 1200.0, 20.0, 300.0, false)'

# After the TW calibration, update the parameter and database then re-run the replay

## vp calibration  used interactive mode.
#root -l 'vpcalib_app.C(0, 0, "vanilla", false, "", "26483,26484,26485,26486,26487,26488", 10.0, "", "./vpcalib_qa", "../../PARAM", -40.0, 40.0)'

# After fitting, produce the pdf using the following command. compareOnly is on. (Or, someone can update the code to produce pdf files at one place. Sorry that I didn't!)
root -l -b -q 'vpcalib_app.C(0, 0, "vanilla", true, "", "26483,26484,26485,26486,26487,26488", 10.0, "", "./vpcalib_qa", "../../PARAM", -40.0, 40.0, true, 0.0, true, true, true)'
#
## lambda calibration can be done in non-interactive mode.
#root -l -b -q 'lambda_calib_app.C(0, "26483,26484,26485,26486,26487,26488", "", true, "", 10.0, 125.0, "", "./lambda_qa", "../../PARAM", "vanilla", true, true, true, false)'


## to only compare lambda between reference (vanilla) and the calibration results, without update, use
root -l -b -q 'lambda_calib_app.C(0, "26483,26484,26485,26486,26487,26488", "", true, "", 10.0, 125.0, "", "./lambda_qa", "../../PARAM", "vanilla", true, true, true, true)'

#After the hodoscope calibration, manually update the param file, for example, for the Vp one,
#for i in {0..99}; do
#   sed -i 's/phodo_Vpcalib/phodo_Vpcalib_26483-26488/g' "general_${i}.param"
#done
# The same thing can be done for TW.
