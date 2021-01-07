"""
Notes
"""

# %%
import os
import subprocess
from gp_step0_dcm2nii import func_sbatch

# %%
work_dir = "/scratch/madlab/nate_vCAT/derivatives"
subj = "sub-005"
sess = "ses-S1"
train = "loc_decon"
test = "Study_BE"

subj_dir = os.path.join(work_dir, subj, sess)
mask = os.path.join(work_dir, "grpAnalysis/Group_Int_Mask.nii.gz")

# %%
h_cmd = f"""
    module load afni-20.2.06
    cd {subj_dir}
    3dsvm -trainvol 3dSVM_{train}+tlrc \
        -trainlabels 3dSVM_{train}_categories.txt \
        -mask {mask} \
        -model 3dSVM_train
"""
h_train = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
train_out = h_train.communicate()[0].decode("utf-8")
# print(out)
f = open(f"{subj_dir}/3dSVM_train_out.txt", "w")
f.write(train_out)
f.close()


# %%
h_cmd = f"""
    module load afni-20.2.06
    cd {subj_dir}
    3dsvm -testvol 3dSVM_{test}+tlrc \
        -model 3dSVM_train+tlrc \
        -testlabels 3dSVM_{test}_categories.txt \
        -predictions 3dSVM_pred > 3dSVM_pred.txt 2>&1
"""
# func_sbatch(h_cmd, 1, 1, 1, f"{subj.split('-')[1]}tst", subj_dir)
h_test = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
test_out = h_test.communicate()

# %%
