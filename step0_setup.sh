#!/bin/bash

parDir=~/Projects/learn_mvpa/vCAT_data

numRuns=4
condList=(hit miss)
taskList=(vCAT)

subjList=(sub006)


# # set up keys
cd $parDir

tr=`3dinfo -tr CleanData_vCAT+tlrc`
echo -e "TR\t$tr" > scan_key.txt
echo "vCAT task" > study_key.txt

for((i=0;i<${#taskList[@]};i++)); do

	pos=$(($i + 1))
	echo -e "task00${pos}\t${taskList[$i]}" > task_key.txt

	keyDir=models/model00$pos
	mkdir -p $keyDir
	print=${keyDir}/condition_key.txt
	> $print

	for j in ${!condList[@]}; do
		c=$(($j + 1))
		echo -e "task00${pos}\tcond00${c}\t${condList[$j]}" >> $print
	done
done


# set up subDir
for i in ${subjList[@]}; do

	# set up dirs
	mkdir -p ${i}/{BOLD,anatomy,model}

	# anat
	3dcopy struct_ns+tlrc ${i}/anatomy/struct_ns.nii.gz

	for j in ${!taskList[@]}; do

		tc=$(($j+1))
		ntim=`3dinfo -ntimes CleanData_${taskList[$j]}+tlrc`
		runLen=$(($ntim / $numRuns))
		start=0
		end=$(($start + $runLen - 1))

		# Split TF into attributes
		for k in ${condList[@]}; do			
			timing_tool.py -timing tf_${taskList[$j]}_${k}.txt -tr $tr \
				-stim_dur 2 -run_len $(echo "$runLen*$tr" | bc) -min_frac 0.3 \
				-timing_to_1D att_${k} -per_run_file
		done

		for((k=1;k<=$numRuns;k++)); do

			## BOLD
			boldDir=${i}/BOLD/task00${tc}_run00${k}
			mkdir -p $boldDir

			# Split CleanData
			3dTcat -prefix tmp_CleanData -tr $tr "CleanData_${taskList[$j]}+tlrc[${start}..${end}]"
			3dcopy tmp_CleanData+tlrc ${boldDir}/bold.nii.gz && rm tmp_CleanData+tlrc.*
			let start+=$runLen
			let end+=$runLen

			# Make Attribute file - this is hardcoded (probably better to do in python)
			att1=(`cat att_hit_r0${k}.1D`)
			att2=(`cat att_miss_r0${k}.1D`)
			print=${boldDir}/attributes.txt
			> $print
			for m in ${!att1[@]}; do
				if [ ${att1[$m]} == 1 ]; then
					echo -e "hit\t0" >> $print
				elif [ ${att2[$m]} == 1 ]; then
					echo -e "miss\t0" >> $print
				else
					echo -e "base\t0" >> $print
				fi
			done

			# model
			modelDir=${i}/model/model00${tc}/onsets/task00${tc}_run00${k}
			mkdir -p $modelDir

			for m in ${!condList[@]}; do
				cc=$(($m+1))
				print=${modelDir}/cond00${cc}.txt
				> $print

				arr=(`sed -n ${cc}p tf_${taskList[$j]}_${condList[$m]}.txt`)
				for n in ${arr[@]}; do
					echo -e "${n}\t2\t1" >> $print
				done
			done
		done
	done
done

