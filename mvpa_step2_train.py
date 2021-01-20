"""
Notes

Updated to train (and then test) on subjects individually.
Group training yielded:
    error = 37.6%
    recall = 25.9% (hit/(hit+miss))
    precision = 92.5% (hit/fa)

    Individual training yields better recall, and lower error.

TODO update to multiple intersection mask x ROI mask
"""

# %%
import os
import sys
from gp_step0_dcm2nii import func_sbatch


# %%
def func_train(subj_str, sess_str, work_dir, train_str):
    """
    Train on data, write model to MVPA_train+tlrc.
    Capture training stats in MVPA_train_acc.txt
    """
    subj_dir = os.path.join(work_dir, subj_str, sess_str)
    if not os.path.exists(os.path.join(subj_dir, "MVPA_train+tlrc.HEAD")):
        h_cmd = f"""
            cd {subj_dir}
            3dsvm -trainvol MVPA_{train_str}+tlrc \
                -trainlabels MVPA_{train_str}_categories.txt \
                -mask mask_epi_anat+tlrc \
                -model MVPA_train \
                > MVPA_train_acc.txt 2>&1
        """
        func_sbatch(h_cmd, 1, 6, 2, f"{subj_str.split('-')[1]}trn", subj_dir)


def main():

    main_subj_str = str(sys.argv[1])
    main_sess_str = str(sys.argv[2])
    main_train_str = str(sys.argv[3])
    main_work_dir = str(sys.argv[4])

    func_train(main_subj_str, main_sess_str, main_work_dir, main_train_str)


if __name__ == "__main__":
    main()

# %%
