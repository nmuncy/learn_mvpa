#!/bin/bash

deriv_dir="/scratch/madlab/nate_vCAT/derivatives"
out_dir=${deriv_dir}/grpAnalysis
if [ ! -d $out_dir ]; then
    mkdir -p $out_dir
fi

print=${out_dir}/Stats_train_acc.txt
echo -e "Subj\tPrec\tRec\tErr" > $print

cd $deriv_dir
for i in sub*; do

    prec=`grep "precision" ${i}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`
    rec=`grep "recall" ${i}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`
    err=`grep "error" ${i}/ses-S1/MVPA_train_acc.txt | awk '{print $6}'`

    echo -e "${i}\t${prec:11}\t${rec:8}\t${err:7}" >> $print
done
