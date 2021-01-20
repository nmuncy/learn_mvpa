import os
import fnmatch
import pandas as pd
import numpy as np


def main():
    """
    This script will combine category and predicted files, split for
        category type, and write out to beh_dir for each test.
    """

    # Set up
    deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    beh_dir = os.path.join(deriv_dir, "grpAnalysis/behAnalysis")
    test_list = ["Study_BE", "Study_CE", "Study_FP"]
    sess_str = "ses-S1"
    cat_dict = {"Face": 1, "Scene": 2}

    # make output dir
    if not os.path.exists(beh_dir):
        os.makedirs(beh_dir)

    # get subjs
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    for subj in subj_list:
        subj_dir = os.path.join(deriv_dir, subj, sess_str)

        for test in test_list:
            if os.path.exists(os.path.join(subj_dir, f"MVPA_pred_{test}.1D")):

                # Combine true, predicted categories
                df_subj = pd.read_table(
                    os.path.join(subj_dir, f"MVPA_{test}_cat_updated.txt")
                )
                df_subj["Pred"] = pd.read_table(
                    os.path.join(subj_dir, f"MVPA_pred_{test}.1D")
                )
                df_subj.columns = ["True", "Pred"]

                # Split into individual cat files, for ROC
                for cat in cat_dict:
                    np.savetxt(
                        os.path.join(beh_dir, f"{subj}_{test}_{cat}.txt"),
                        df_subj[
                            df_subj["True"].astype(str).str.contains(f"{cat_dict[cat]}")
                        ],
                        fmt="%s",
                        delimiter=",",
                        header="True,Pred",
                        comments="",
                    )


if __name__ == "__main__":
    main()
