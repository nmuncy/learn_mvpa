"""
Notes

Updated to use individual subject training models.
"""


import os
import numpy as np
import pandas as pd
import sys
import json
from gp_step0_dcm2nii import func_sbatch


def func_test(subj, test_list, subj_dir):

    """
    Generates a test file consisting of only relevant volumes.
        This does not change the output, but I like it and it will
        be useful when splitting the phase into runs.

    Then, the classifier attempts to predict those volumes. The output files
        are categorical (class), and afni writes stats to stderr.
    """

    subj_num = subj.split("-")[-1]

    # support multiple decons
    for test in test_list:

        # Get only test volumes
        if not os.path.exists(os.path.join(subj_dir, f"MVPA_{test}_test+tlrc.HEAD")):

            # make new category file for only relevant volumes
            #   will be used for ROC
            df = pd.read_csv(
                os.path.join(subj_dir, f"MVPA_{test}_categories.txt"),
                sep=" ",
                header=None,
            )
            df.columns = ["cat"]
            df_new = df[df["cat"] != 9999]
            df_out = os.path.join(subj_dir, f"MVPA_{test}_cat_updated.txt")
            np.savetxt(df_out, df_new["cat"].values, fmt="%s", delimiter=" ")

            # make file of only relevant volumes
            df_list = df_new.index.values.tolist()
            h_cmd = f"""
                cd {subj_dir}
                3dTcat -prefix MVPA_{test}_test MVPA_{test}+tlrc[{",".join(str(i) for i in df_list)}]
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}cat", subj_dir)

        # Test
        if not os.path.exists(os.path.join(subj_dir, f"MVPA_pred_{test}.1D")):
            h_cmd = f"""
                cd {subj_dir}
                3dsvm \
                    -testvol MVPA_{test}_test+tlrc \
                    -model MVPA_train+tlrc \
                    -testlabels MVPA_{test}_cat_updated.txt \
                    -predictions MVPA_pred_{test} \
                    -classout \
                    2> MVPA_pred_{test}_acc.txt
            """
            out_str = f"{subj_num}{test.split('_')[1]}"
            func_sbatch(h_cmd, 1, 8, 4, out_str, subj_dir)


def main():
    main_subj = str(sys.argv[1])
    main_subj_dir = str(sys.argv[2])
    with open(os.path.join(main_subj_dir, "test_list.json")) as json_file:
        main_test_list = json.load(json_file)

    func_test(main_subj, main_test_list, main_subj_dir)


if __name__ == "__main__":
    main()
