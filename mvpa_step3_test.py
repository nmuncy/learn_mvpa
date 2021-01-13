"""
Notes
"""


import os
import numpy as np
import pandas as pd
import sys
import json
import subprocess
import time

# from gp_step0_dcm2nii import func_sbatch


def func_sbatch(command, wall_hours, mem_gig, num_proc, h_str, work_dir):
    #
    full_name = f"{work_dir}/sbatch_writeOut_{h_str}"
    sbatch_job = f"""
        sbatch \
        -J {h_str} -t {wall_hours}:00:00 --mem={mem_gig}000 --ntasks-per-node={num_proc} \
        -p IB_40C_512G -o {full_name}.out -e {full_name}.err \
        --account iacc_madlab --qos pq_madlab \
        --wrap="module load afni-20.2.06 \n {command}"
    """
    sbatch_response = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_response.communicate()[0]
    print(job_id, h_str, sbatch_job)

    while_count = 0
    status = False
    while not status:

        check_cmd = "squeue -u $(whoami)"
        sq_check = subprocess.Popen(check_cmd, shell=True, stdout=subprocess.PIPE)
        out_lines = sq_check.communicate()[0]
        b_decode = out_lines.decode("utf-8")

        if h_str not in b_decode:
            status = True

        if not status:
            while_count += 1
            print(f"Wait count for sbatch job {h_str}: ", while_count)
            time.sleep(3)
    print(f'Sbatch job "{h_str}" finished')


def func_test(subj, test_list, subj_dir, group_dir):

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
        if not os.path.exists(os.path.join(subj_dir, f"3dSVM_{test}_test+tlrc.HEAD")):
            df = pd.read_csv(
                os.path.join(subj_dir, f"3dSVM_{test}_categories.txt"),
                sep=" ",
                header=None,
            )
            df.columns = ["cat"]

            df_new = df[df["cat"] != 9999]
            df_out = os.path.join(subj_dir, f"3dSVM_{test}_cat_updated.txt")
            np.savetxt(df_out, df_new["cat"].values, fmt="%s", delimiter=" ")

            df_list = df_new.index.values.tolist()
            h_cmd = f"""
                    cd {subj_dir}
                    3dTcat -prefix 3dSVM_{test}_test 3dSVM_{test}+tlrc[{",".join(str(i) for i in df_list)}]
                """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}cat", subj_dir)

        # Test
        if not os.path.exists(os.path.join(subj_dir, f"3dSVM_pred_{test}.1D")):
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
            func_sbatch(h_cmd, 1, 4, 4, f"{subj_num}tst", subj_dir)


def main():
    main_subj = str(sys.argv[1])
    main_subj_dir = str(sys.argv[2])
    main_group_dir = str(sys.argv[3])
    with open(os.path.join(main_subj_dir, "test_list.json")) as json_file:
        main_test_list = json.load(json_file)

    func_test(main_subj, main_test_list, main_subj_dir, main_group_dir)


if __name__ == "__main__":
    main()
