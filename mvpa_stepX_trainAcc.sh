#!/bin/bash

parDir=/scratch/madlab/nate_vCAT
derivDir=${parDir}/derivatives
outDir=${parDir}/analyses

print=${outDir}/Stats_train_acc.txt
echo "Note: recall = hit/(hit+miss); precision = hit/fa" > $print
echo >> $print
echo -e "Subj \t Prec \t Rec \t Err" >> $print

cd $derivDir
for subj in sub*; do

    prec=`grep "precision" ${subj}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`
    rec=`grep "recall" ${subj}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`
    err=`grep "error" ${subj}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`

    echo -e "${subj} \t ${prec:11} \t ${rec:8} \t ${err:7}" >> $print
done
