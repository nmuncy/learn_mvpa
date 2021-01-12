"""
Notes
"""

# %%
import os
import sys
import fnmatch
from gp_step0_dcm2nii import func_sbatch


# %%
def func_train(subj_list, sess_str, train_str, work_dir, group_dir, mask):

    # Combine all train epi and category files
    train_list = []
    txt_list = []
    for subj in subj_list:
        h_file = os.path.join(work_dir, subj, sess_str, f"3dSVM_{train_str}+tlrc")
        h_txt = os.path.join(
            work_dir, subj, sess_str, f"3dSVM_{train_str}_categories.txt"
        )
        if os.path.exists(f"{h_file}.HEAD"):
            train_list.append(h_file)
            txt_list.append(h_txt)

    if not os.path.exists(os.path.join(group_dir, "vCAT_train+tlrc.HEAD")):
        h_cmd = f"""
            cd {group_dir}
            3dTcat -prefix vCAT_train {" ".join(train_list)}
        """
        func_sbatch(h_cmd, 1, 1, 1, "MVPAcomb", group_dir)

    with open(os.path.join(group_dir, "vCAT_train_cat.txt"), "w") as outfile:
        for fname in txt_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    # Train
    if not os.path.exists(os.path.join(group_dir, "3dSVM_train+tlrc.HEAD")):
        h_cmd = f"""
            cd {group_dir}
            3dsvm -trainvol vCAT_train+tlrc \
                -trainlabels vCAT_train_cat.txt \
                -mask {mask} \
                -model 3dSVM_train \
                > 3dSVM_train_acc.txt 2>&1
        """
        func_sbatch(h_cmd, 4, 6, 4, "MVPAtrn", group_dir)


def main():

    main_sess_str = str(sys.argv[1])
    main_train_str = str(sys.argv[2])
    main_work_dir = str(sys.argv[3])
    main_group_dir = str(sys.argv[4])
    main_mask = str(sys.argv[5])

    main_subj_list = [
        x for x in os.listdir(main_work_dir) if fnmatch.fnmatch(x, "sub-*")
    ]
    main_subj_list.sort()

    print(
        main_subj_list,
        main_sess_str,
        main_train_str,
        main_work_dir,
        main_group_dir,
        main_mask,
    )

    func_train(
        main_subj_list,
        main_sess_str,
        main_train_str,
        main_work_dir,
        main_group_dir,
        main_mask,
    )


if __name__ == "__main__":
    main()

# %%
