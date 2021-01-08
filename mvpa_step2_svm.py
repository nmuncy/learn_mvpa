"""
Notes
"""

# %%
import os
import subprocess
import pandas as pd
import numpy as np

# from gp_step0_dcm2nii import func_sbatch

# %%
# Set up
work_dir = "/scratch/madlab/nate_vCAT/derivatives"
subj = "sub-005"
sess = "ses-S1"
train = "loc_decon"
test_list = ["Study_BE", "Study_BP", "Study_CP", "Study_FP"]

subj_dir = os.path.join(work_dir, subj, sess)
mask = os.path.join(work_dir, "grpAnalysis/Group_Int_Mask.nii.gz")

# %%
# Train
if not os.path.exists(os.path.join(subj_dir, "3dSVM_train+tlrc.HEAD")):
    h_cmd = f"""
        module load afni-20.2.06
        cd {subj_dir}
        3dsvm -trainvol 3dSVM_{train}+tlrc \
            -trainlabels 3dSVM_{train}_categories.txt \
            -mask {mask} \
            -model 3dSVM_train \
            > 3dSVM_train_acc.txt 2>&1
    """
    h_train = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    train_out = h_train.communicate()[0].decode("utf-8")
    # f = open(f"{subj_dir}/3dSVM_train_out.txt", "w")
    # f.write(train_out)
    # f.close()


# %%

for test in test_list:
    # test = "Study_BE"

    # Get only test sub-bricks
    #   doesn't affect output, but I like it
    #   TODO split into individual runs

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
                -model 3dSVM_train+tlrc \
                -testlabels 3dSVM_{test}_cat_updated.txt \
                -predictions 3dSVM_pred_{test} \
                -classout \
                2> 3dSVM_pred_{test}_acc.txt
        """
        # func_sbatch(h_cmd, 1, 1, 1, f"{subj.split('-')[1]}tst", subj_dir)
        h_test = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        test_out = h_test.communicate()
        print(test, test_out)

# %%
