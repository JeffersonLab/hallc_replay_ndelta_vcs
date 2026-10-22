# Each app supports interactive/ noninteractive mode, but only noninteractive mode was used for TW and lambda part
# root -l -b -q 'timewalk_calib_app.C(26138, 0, "vanilla", true, "0", "", 10.0)'
# root -l -b -q 'timewalk_calib_app.C(26927, 0, "vanilla", true, "1", "", 10.0)'
# compareOnly is now true, so it will not update the TW fit. compareOnly=kFALSE can be used as an argument, or one can update timewalk_calib_app.C

## vp calibration  used interactive mode.
root -l 'vpcalib_app.C(26138, 0, "vanilla", false, "0", "", 10.0)'
root -l 'vpcalib_app.C(26927, 0, "vanilla", false, "1", "", 10.0, "", "./vpcalib_qa", "../../PARAM", false, 0.0, false, false)'

# After fitting, produce the pdf using the following command. compareOnly is on. (Or, someone can update the code to produce pdf files at one place. Sorry that I didn't!)
root -l -b -q 'vpcalib_app.C(26138, 0, "vanilla", true, "0", "", 10.0, "", "./vpcalib_qa", "../../PARAM", true, 0.0, false, true, true, true)'
root -l -b -q 'vpcalib_app.C(26927, 0, "vanilla", true, "1", "", 10.0, "", "./vpcalib_qa", "../../PARAM", false, 0.0, false, false, true, true)'

## lambda calibration
#root -l -b -q 'lambda_calib_app.C(26138, "", "0", true, "0", 10.0, 100.0, "", "./lambda_qa", "../../PARAM", "vanilla", true, true, true, false)'
#root -l -b -q 'lambda_calib_app.C(26927, "", "1", true, "1", 10.0, 100.0, "", "./lambda_qa", "../../PARAM", "vanilla", false, false, true, false)'

# for comparison only, without updating the parameter file, use the compareOnly flag
root -l -b -q 'lambda_calib_app.C(26138, "", "0", true, "0", 10.0, 100.0, "", "./lambda_qa", "../../PARAM", "vanilla", true, true, true, true)'
root -l -b -q 'lambda_calib_app.C(26927, "", "1", true, "1", 10.0, 100.0, "", "./lambda_qa", "../../PARAM", "vanilla", false, false, true, true)'

# After the calibration, manually update the TW and VP, ex)
# for i in {0..81}; do   sed -i '' 's/hhodo_Vpcalib/hhodo_Vpcalib_0/g' "general_${i}.param"; done
# for i in {82..99}; do   sed -i '' 's/hhodo_Vpcalib/hhodo_Vpcalib_1/g' "general_${i}.param"; done
# The same thing can be done for TW.
