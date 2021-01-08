"""
Notes
"""

# %%
import os
import subprocess
import pandas as pd
import numpy as np
import fnmatch
from gp_step0_dcm2nii import func_sbatch


# %%
# Set up
work_dir = "/scratch/madlab/nate_vCAT/derivatives"
group_dir = os.path.join(work_dir, "grpAnalysis")
sess = "ses-S1"
train = "loc_decon"
test_list = ["Study_BE", "Study_BP", "Study_CP", "Study_FP"]
mask = os.path.join(group_dir, "Group_Int_Mask.nii.gz")


# %%
# Combine all train data
subj_list = [x for x in os.listdir(work_dir) if fnmatch.fnmatch(x, "sub-*")]
subj_list.sort()

train_list = []
txt_list = []
for subj in subj_list:
    h_file = os.path.join(work_dir, subj, sess, f"3dSVM_{train}+tlrc")
    h_txt = os.path.join(work_dir, subj, sess, f"3dSVM_{train}_categories.txt")
    if os.path.exists(f"{h_file}.HEAD"):
        train_list.append(h_file)
        txt_list.append(h_txt)

if not os.path.exists(os.path.join(group_dir, "vCAT_train+tlrc.HEAD")):
    h_cmd = f"""
        module load afni-20.2.06
        cd {group_dir}
        3dTcat -prefix vCAT_train {" ".join(train_list)}
    """
    func_sbatch(h_cmd, 1, 1, 1, "MVPAcomb", group_dir)

with open(os.path.join(group_dir, "vCAT_train_cat.txt"), "w") as outfile:
    for fname in txt_list:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)


# %%
# Train
if not os.path.exists(os.path.join(group_dir, "3dSVM_train+tlrc.HEAD")):
    h_cmd = f"""
        module load afni-20.2.06
        cd {group_dir}
        3dsvm -trainvol vCAT_train+tlrc \
            -trainlabels vCAT_train_cat.txt \
            -mask {mask} \
            -model 3dSVM_train \
            > 3dSVM_train_acc.txt 2>&1
    """
    func_sbatch(h_cmd, 1, 4, 4, "MVPAtrn", group_dir)

    # h_train = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    # train_out = h_train.communicate()[0].decode("utf-8")
    # # f = open(f"{subj_dir}/3dSVM_train_out.txt", "w")
    # f.write(train_out)
    # f.close()


# %%
subj = "sub-005"
subj_dir = os.path.join(work_dir, subj, sess)


# for test in test_list:
test = "Study_BE"

# Get only test sub-bricks
#   doesn't affect output, but I like it
if not os.path.exists(os.path.join(subj_dir, f"3dSVM_{test}_test+tlrc.HEAD")):
    df = pd.read_csv(
        os.path.join(subj_dir, f"3dSVM_{test}_categories.txt"), sep=" ", header=None
    )
    df.columns = ["cat"]

    df_new = df[df["cat"] != 9999]
    df_out = os.path.join(subj_dir, f"3dSVM_{test}_cat_updated.txt")
    np.savetxt(df_out, df_new["cat"].values, fmt="%s", delimiter=" ")

    df_list = df_new.index.values.tolist()
    h_cmd = f"""
        module load afni-20.2.06
        cd {subj_dir}
        3dTcat -prefix 3dSVM_{test}_test 3dSVM_{test}+tlrc[{",".join(str(i) for i in df_list)}]
    """
    h_cat = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_cat.communicate()

# Test
if not os.path.exists(os.path.join(subj_dir, f"3dSVM_pred_{test}_acc.txt")):
    h_cmd = f"""
        module load afni-20.2.06
        cd {subj_dir}
        3dsvm -testvol 3dSVM_{test}_test+tlrc \
            -model {group_dir}/3dSVM_train+tlrc \
            -testlabels 3dSVM_{test}_cat_updated.txt \
            -predictions 3dSVM_pred_{test} \
            -classout \
            2> 3dSVM_pred_{test}_acc.txt
    """
    h_test = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    test_out = h_test.communicate()
    print(test, test_out)

# %%
