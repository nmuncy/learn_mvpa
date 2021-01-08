"""
Notes
"""

# %%
import os
import subprocess

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
# Train, capture stdout/err which contains stats
if not os.path.exists(os.path.join(subj_dir, "3dSVM_train+tlrc.HEAD")):
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
    f = open(f"{subj_dir}/3dSVM_train_out.txt", "w")
    f.write(train_out)
    f.close()


# %%
# Test, capture stderr for accuracy report
for test in test_list:
    if not os.path.exists(os.path.join(subj_dir, f"3dSVM_pred_{test}_acc.txt")):
        h_cmd = f"""
            module load afni-20.2.06
            cd {subj_dir}
            3dsvm -testvol 3dSVM_{test}+tlrc \
                -model 3dSVM_train+tlrc \
                -testlabels 3dSVM_{test}_categories.txt \
                -predictions 3dSVM_pred_{test} \
                -classout \
                2> 3dSVM_pred_{test}_acc.txt
        """
        # func_sbatch(h_cmd, 1, 1, 1, f"{subj.split('-')[1]}tst", subj_dir)
        h_test = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
        test_out = h_test.communicate()
        print(test, test_out)

# %%
