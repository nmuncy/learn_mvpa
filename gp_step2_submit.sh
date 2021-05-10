#!/bin/bash

module load r-3.5.1-gcc-8.2.0-djzshna

behavDir=/home/data/madlab/Mattfeld_vCAT/behav
derivDir=/scratch/madlab/nate_vCAT/derivatives
codeDir=~/compute/learn_mvpa
phaseArr=("loc" "task")
runArr=(2 4)

# # For testing
# behavDir=/Users/nmuncy/Projects/learn_mvpa/vCAT_data
# derivDir=/Users/nmuncy/Projects/afni_python
# codeDir=$derivDir
# numRuns=2

cd $behavDir
for i in vCAT*; do

	subj=sub-${i#*_}
	dataDir=${behavDir}/$i
	outDir=${derivDir}/${subj}/ses-S1

	if [ -d $outDir ]; then
		for j in ${!phaseArr[@]}; do
			Rscript ${codeDir}/gp_step2_timingFiles.R $dataDir $outDir $i ${runArr[$j]} ${phaseArr[$j]}
		done
	fi

	# write key for task tfs
	print=${outDir}/key_Study.txt
	cat > $print << EOF
Key for timing files tf_Study_*. 
Foo: Letter, Description
	Foo = timing file (tf_study) identifier and decon beh
	Letter = corresponding letter in experiment diagram

Bfe: E, Baseline face event
Bfp: D, Baseline face event and preceding fixed trial
Bse: E, Baseline scene event
Bsp: D, Baseline scene event and preceding fixed trial

Ffpc: F, Fixed face preceding correct conditional trial
Ffpi: F, Fixed face preceding incorrect conditional trial
Fspc: F, Fixed scene preceding correct conditional trial
Fspi: F, Fixed scene preceding incorrect conditional trial

cfec: A, Correct conditional face event
cfei: A, Incorrect conditional face event
csec: A, Correct conditional scene event
csei: A, Incorrect conditional scene event

cfpc: B, Correct conditional face and preceding fixed trial
cfpi: B, Incorrect conditional face and preceding fixed trial
cspc: B, Correct conditional scene and preceding fixed trial
cspi: B, Incorrect conditional scene and preceding fixed trial
EOF

done
