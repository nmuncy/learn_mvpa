#!/bin/bash

parDir=/scratch/madlab/nate_vCAT
workDir=${parDir}/derivatives
checkFile=${parDir}/analyses/FINAL_loc_face-scene+tlrc.HEAD

if [ -f $checkFile ]; then
    cd $workDir
    rm sub-*/ses-S1/*tent*
fi