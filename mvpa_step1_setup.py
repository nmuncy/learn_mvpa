"""
Notes

Written in Python 3.8

Assumes gp_step1-3 from afni_python have already been done
"""

# %%
import subprocess
import fnmatch
import os
import json
import time
import pandas as pd
import numpy as np
from shutil import copyfile
from argparse import ArgumentParser
from gp_step0_dcm2nii import func_sbatch


def func_detrend(tmp_dir, decon_str, tmp_list, tmp_tr):

    # get relevant column numbers
    #   read design matrix, find column labels,
    #   then look for column label in task_dict
    #   to get a brick_list of only effects of interest
    brick_list = []
    with open(os.path.join(tmp_dir, f"X.{decon_str}.xmat.1D")) as f:
        h_file = f.readlines()
        for line in h_file:
            if line.__contains__("ColumnLabels"):
                col_list = line.split('"')[1].replace(" ", "").split(";")
                for i, j in enumerate(col_list):
                    if j.split("#")[0] in tmp_list:
                        brick_list.append(f"{str(i)}")

    # determine correct number of sub-bricks
    #   decon adds an extra for model
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -nv {tmp_dir}/{decon_str}_cbucket_REML+tlrc"
    h_len = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    len_wrong = h_len.communicate()[0].decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # make detrended data
    #   only contains tent data, extra sub-bricks
    #   are left behind
    h_cmd = f"""
        cd {tmp_dir}
        3dTcat -prefix tmp_{decon_str}_cbucket -tr {tmp_tr} "{decon_str}_cbucket_REML+tlrc[0..{len_right}]"
        3dSynthesize -prefix MVPA_{decon_str}_all -matrix X.{decon_str}.xmat.1D \
                -cbucket tmp_{decon_str}_cbucket+tlrc -select {" ".join(brick_list)} -cenfill nbhr
    """
    func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}all", tmp_dir)


def func_lenRun(tmp_dir, decon_str, phase):

    # determine number of volumes
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -ntimes {subj_dir}/MVPA_{decon_str}_all+tlrc"
    h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    num_nvol = int(h_nvol.communicate()[0].decode("utf-8").strip())

    # determine number of runs
    num_runs = len(
        [
            x
            for x in os.listdir(subj_dir)
            if fnmatch.fnmatch(x, f"*{phase}*scale+tlrc.HEAD")
        ]
    )

    # determine length of e/run
    len_run = int(num_nvol / num_runs)
    tmp_dict = {"NVol": num_nvol, "NRun": num_runs, "LRun": len_run}
    return tmp_dict


def func_split(tmp_dir, phase, decon_str, task_num):

    # split mvpa file into individual runs
    h_dict = func_lenRun(tmp_dir, decon_str, phase)
    len_run = h_dict["LRun"]

    beg_vol = 0
    end_vol = len_run - 1
    for run in range(1, h_dict["NRun"] + 1):

        # make dir for e/run
        bold_dir = os.path.join(mvpa_dir, f"BOLD/task00{task_num}_run00{run}")
        if not os.path.exists(bold_dir):
            os.makedirs(bold_dir)

        # split
        if not os.path.exists(os.path.join(bold_dir, "bold.nii.gz")):
            h_cmd = f"""
                cd {tmp_dir}
                3dTcat -prefix tmp_run-{run}_{decon_str}_MVPA -tr {len_tr} "MVPA_{decon_str}_all+tlrc[{beg_vol}..{end_vol}]"
                3dcopy tmp_run-{run}_{decon_str}_MVPA+tlrc {bold_dir}/bold.nii.gz
            """
            func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}spl", tmp_dir)
        beg_vol += len_run
        end_vol += len_run


# %%
# def func_job(subj, subj_dir, len_tr, task_dict, der_dir):
"""
Step 1: Detrend

Extract TENT sub-bricks associated with behaviors
    into a single file
"""
# For testing
subj = "sub-005"
subj_dir = "/scratch/madlab/nate_vCAT/derivatives/sub-005/ses-S1"
len_tr = 1.76
task_dict = {
    "loc": ["face", "scene", "num"],
    "Study": {"BE": ["Bfe", "Bse"], "FP": ["Ffpc", "Ffpi", "Fspc", "Fspi"]},
}
beh_dur = 1
der_dir = "/scratch/madlab/nate_vCAT/derivatives"


# %%
# Work
subj_num = subj.split("-")[1]
mvpa_dir = os.path.join(der_dir, f"mvpa/sub{subj_num}")

for phase in task_dict:
    if type(task_dict[phase]) == list:
        decon_str = f"{phase}_decon"
        if not os.path.exists(
            os.path.join(subj_dir, f"MVPA_{decon_str}_all+tlrc.HEAD")
        ):
            func_detrend(subj_dir, decon_str, task_dict[phase], len_tr)
    elif type(task_dict[phase]) == dict:
        for decon in task_dict[phase]:
            decon_str = f"{phase}_{decon}"
            if not os.path.exists(
                os.path.join(subj_dir, f"MVPA_{decon_str}_all+tlrc.HEAD")
            ):
                func_detrend(subj_dir, decon_str, task_dict[phase][decon], len_tr)


# %%
"""
Step 2: Organize MRI data

1) Pymvpa expects certain files in certain locations
2) Split MVPA afni file into individual runs
"""
# pymvpa dirs
py_dirs = ["BOLD", "anatomy", "model", "masks"]
for i in py_dirs:
    if not os.path.exists(os.path.join(mvpa_dir, i)):
        os.makedirs(os.path.join(mvpa_dir, i))

# anat
if not os.path.exists(os.path.join(mvpa_dir, "anatomy/struct_ns.nii.gz")):
    h_cmd = f"module load afni-20.2.06 \n 3dcopy {subj_dir}/struct_ns+tlrc {mvpa_dir}/anatomy/struct_ns.nii.gz"
    h_job = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    print(h_job.communicate())

# masks
mask_dir = os.path.join(mvpa_dir, "masks", "orig")
if not os.path.exists(mask_dir):
    os.makedirs(mask_dir)

if not os.path.exists(os.path.join(mask_dir, "Group_Int_Mask.nii.gz")):
    copyfile(
        os.path.join(der_dir, "grpAnalysis/Group_Int_Mask.nii.gz"),
        os.path.join(mask_dir, "Group_Int_Mask.nii.gz"),
    )

# BOLD - split into runs
task_count = 1
for phase in task_dict:
    if type(task_dict[phase]) == list:
        h_str = f"{phase}_decon"
        func_split(subj_dir, phase, h_str, task_count)
        task_count += 1
    elif type(task_dict[phase]) == dict:
        for decon in task_dict[phase]:
            h_str = f"{phase}_{decon}"
            func_split(subj_dir, phase, h_str, task_count)
            task_count += 1


# %%
"""
Step 3: Organize timing files

1) Based on timing files (e.g. tf_loc_face.txt) from
    gp_step2_timingFiles.R. Will generate needed attribute files.

2) attributes.txt - one attribute per volume
    - written to sub*/BOLD/task*

    - Issue from behavior duration, TR mismatch
    resulting in holes.

    - Fix - fill NaN with preceding attribute value
    since study is a block design.

3) cond00?.txt - each condition/attribut type has own onset time
    - ref condtion_key.txt file for definitions
    - written to sub*/model/model00?/onsets/task*
    - format = onset, block duration, 1
        e.g. 157.5 22.5 1
        duration is in volume time
"""

# %%


def func_timing(beh_list, len_tr, tmp_dir, phase, decon_str):

    # determine phase length in seconds
    h_dict = func_lenRun(tmp_dir, decon_str, phase)
    len_sec = h_dict["LRun"] * len_tr

    # split timing files into 1D
    for beh in beh_list:
        if not os.path.exists(os.path.join(subj_dir, f"tmp_tf_{phase}_{beh}_r01.1D")):

            # determine behavior duration, get 1st number in column 0
            data_raw = pd.read_csv(
                os.path.join(tmp_dir, f"dur_{phase}_{beh}.txt"), header=None
            )
            data_clean = data_raw[0].str.split("\t", expand=True)
            for value in data_clean[0]:
                if "*" not in value:
                    beh_dur = value
                    break

            # make 1D file per run
            h_cmd = f"""
                module load afni-20.2.06
                timing_tool.py -timing {subj_dir}/tf_{phase}_{beh}.txt \
                    -tr {len_tr} -stim_dur {beh_dur} -run_len {len_sec} -timing_to_1D \
                    {subj_dir}/tmp_tf_{phase}_{beh} -per_run_file
            """
            h_spl = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            print(h_spl.communicate())


for count, phase in enumerate(task_dict):

    if type(task_dict[phase]) == list:
        func_timing(task_dict[phase], len_tr, subj_dir, phase, f"{phase}_decon")

    # %%
    # make attribute files for e/run
    for run in range(1, num_runs + 1):

        # determine split tf files
        tf_list = [
            x
            for x in os.listdir(subj_dir)
            if fnmatch.fnmatch(x, f"tmp_tf_{phase}*_r0{run}.1D")
        ]
        tf_list.sort()

        # make attribute df - combine columns
        df_att = pd.concat(
            [
                pd.read_csv(os.path.join(subj_dir, item), names=[item.split("_")[3]])
                for item in tf_list
            ],
            axis=1,
        )
        column_list = list(df_att)

        # new attribute (att) column, fill with column name (cond)
        for cond in cond_list:
            df_att.loc[df_att[cond] == 1, "att"] = cond

        # clean up df
        for i, j in enumerate(df_att["att"]):

            # Fill holes with preceding att value
            #   skip first, last few volumes
            if i > 1 and i < (len(df_att) - 2) and pd.isnull(df_att.loc[i, "att"]):
                if not pd.isnull(df_att.loc[(i - 1), "att"]):
                    df_att.at[i, "att"] = df_att.loc[(i - 1), "att"]
                else:
                    df_att.at[i, "att"] = df_att.loc[(i - 2), "att"]

            # remove volumes that had multiple behaviors
            if df_att[column_list].iloc[i].sum() > 1:
                df_att.at[i, "att"] = "base"

        # fill NaN with rest, add col of 0s for some reason
        df_att = df_att.replace(np.nan, "base", regex=True)
        df_att["zero"] = 0

        # write
        col_select = ["att", "zero"]
        h_out = os.path.join(
            mvpa_dir, "BOLD", f"task00{count+1}_run00{run}", "attributes.txt"
        )
        np.savetxt(h_out, df_att[col_select].values, fmt="%s", delimiter=" ")

        # make onset times for each task/run
        model_dir = os.path.join(
            mvpa_dir,
            "model",
            f"model00{count+1}",
            "onsets",
            f"task00{count+1}_run00{run}",
        )
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        for cc, cond in enumerate(cond_list):

            # determine start and end volume of each cond block
            ons_dict = {"Start": [], "End": []}
            for i, j in enumerate(df_att["att"]):
                if j == cond and not df_att.loc[(i - 1), "att"] == cond:
                    ons_dict["Start"].append(i)
                elif j == cond and not df_att.loc[(i + 1), "att"] == cond:
                    ons_dict["End"].append(i)

            # determine start, duration of each block in volume time
            df_ons = pd.DataFrame(ons_dict)
            df_ons["ons"] = round((df_ons["Start"] * len_tr), 1)
            df_ons["dur"] = round(((df_ons["End"] - df_ons["Start"]) * len_tr), 1)
            df_ons["one"] = 1

            # write
            col_select = ["ons", "dur", "one"]
            h_out = os.path.join(model_dir, f"cond00{cc+1}.txt")
            np.savetxt(h_out, df_ons[col_select].values, fmt="%s", delimiter=" ")


# %%
# # receive arguments
# def func_argparser():
#     parser = ArgumentParser("Receive Bash args from wrapper")
#     parser.add_argument("h_sub", help="Subject ID")
#     parser.add_argument("h_dir", help="Subject Directory")
#     parser.add_argument("h_trl", help="TR Length")
#     parser.add_argument("h_der", help="Derivatives Directory")
#     return parser


# # %%
# def main():

#     args = func_argparser().parse_args()
#     with open(os.path.join(args.h_der, "mvpa/task_dict.json")) as json_file:
#         h_task_dict = json.load(json_file)

#     # print(args.h_sub, args.h_dir, args.h_trl, h_task_dict, args.h_der)
#     func_job(
#         args.h_sub,
#         args.h_dir,
#         float(args.h_trl),
#         h_task_dict,
#         args.h_der,
#     )


# if __name__ == "__main__":
#     main()

# %%
