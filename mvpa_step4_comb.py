# %%
import os
import fnmatch
import pandas as pd
import numpy as np


# %%
def main():
    """
    This script will combine category and predicted files, split for
        category type, and write out to beh_dir for each test.
    """

    # Set up
    par_dir = "/scratch/madlab/nate_vCAT"
    deriv_dir = os.path.join(par_dir, "derivatives")
    beh_dir = os.path.join(par_dir, "analyses/mvpa_pred")
    test_list = ["Study_single"]
    sess_str = "ses-S1"

    # you can determine cat_dict from reviewing the
    #   MVPA_*_matrix.txt files
    cat_dict = {
        "Face": 1,
        "Scene": 2,
    }

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
                    os.path.join(subj_dir, f"MVPA_{test}_cat_updated.txt"),
                    header=None,
                )
                df_subj["Pred"] = pd.read_table(
                    os.path.join(subj_dir, f"MVPA_pred_{test}.1D"),
                    header=None,
                )
                df_subj.columns = ["True", "Pred"]

                # save all
                df_subj.to_csv(
                    os.path.join(beh_dir, f"{subj}_{test}_all.txt"),
                    index=False,
                    encoding="utf-8",
                )

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
