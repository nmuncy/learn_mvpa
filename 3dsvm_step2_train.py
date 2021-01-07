"""
Notes
"""

# %%
import os
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
if not os.path.exists(os.path.join(subj_dir, "3dSVM_train+tlrc.HEAD")):
    h_cmd = f"""
        cd {subj_dir}
        3dsvm -trainvol 3dSVM_{train}+tlrc \
            -trainlabels 3dSVM_{train}_categories.txt \
            -mask {mask} \
            -model 3dSVM_train
    """
    func_sbatch(h_cmd, 1, 1, 1, f"{subj.split('-')[1]}trn", subj_dir)


# %%
if not os.path.exists(os.path.join(subj_dir, "3dSVM_test_overall_DAG.1D")):
    h_cmd = f"""
        cd {subj_dir}
        3dsvm -testvol 3dSVM_{test}+tlrc \
            -model 3dSVM_train+tlrc \
            -testlabels 3dSVM_{test}_categories.txt \
            -predictions 3dSVM_test
    """
    func_sbatch(h_cmd, 1, 1, 1, f"{subj.split('-')[1]}tst", subj_dir)


# %%
