"""
Notes
"""

# %%
import os
import fnmatch
import subprocess
import pandas as pd
import numpy as np
from gp_step0_dcm2nii import func_sbatch


# %%
def func_lenRun(lrn_subj_dir, lrn_phase):

    # determine number of volumes
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -ntimes {lrn_subj_dir}/run-1_{lrn_phase}_scale+tlrc"
    h_nvol = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    num_nvol = int(h_nvol.communicate()[0].decode("utf-8").strip())

    # determine TR len
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -tr {lrn_subj_dir}/run-1_{lrn_phase}_scale+tlrc"
    h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    len_tr = float(h_tr.communicate()[0].decode("utf-8").strip())

    # determine number of runs
    num_runs = len(
        [
            x
            for x in os.listdir(lrn_subj_dir)
            if fnmatch.fnmatch(x, f"*{lrn_phase}*scale+tlrc.HEAD")
        ]
    )

    # determine length of e/run
    len_run = float(num_nvol * len_tr)

    tmp_dict = {
        "LenTR": len_tr,
        "NumVol": num_nvol,
        "NumRun": num_runs,
        "LenRun": len_run,
    }
    return tmp_dict


# %%
def func_timing(
    tim_beh_list,
    tim_hdr_dict,
    tim_subj_dir,
    tim_phase,
    tim_dcn_str,
):
    """
    Make attributes files
    """
    # # For Testing
    # tim_beh_list = task_dict["loc"]
    # tim_hdr_dict = hdr_dict
    # tim_subj_dir = subj_dir
    # tim_phase = "loc"

    # convert timing files to 1D in volume time
    for beh in tim_beh_list:
        if not os.path.exists(
            os.path.join(tim_subj_dir, f"tmp_tf_{tim_phase}_{beh}.1D")
        ):

            # determine behavior duration, get 1st number in column 0
            data_raw = pd.read_csv(
                os.path.join(tim_subj_dir, f"dur_{tim_phase}_{beh}.txt"), header=None
            )
            data_clean = data_raw[0].str.split("\t", expand=True)
            for value in data_clean[0]:
                if "*" not in value:
                    beh_dur = value
                    break

            # make 1D file per run
            h_cmd = f"""
                module load afni-20.2.06
                timing_tool.py -timing {tim_subj_dir}/tf_{tim_phase}_{beh}.txt \
                    -tr {tim_hdr_dict["LenTR"]} -stim_dur {beh_dur} \
                    -run_len {tim_hdr_dict["LenRun"]} -timing_to_1D \
                    {tim_subj_dir}/tmp_tf_{tim_phase}_{beh}.1D
            """
            h_spl = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
            print(h_spl.communicate())

    # make attribute df - combine columns
    tf_list = []
    for beh in tim_beh_list:
        h_tf = f"tmp_tf_{tim_phase}_{beh}.1D"
        if os.path.exists(os.path.join(tim_subj_dir, h_tf)):
            tf_list.append(h_tf)

    df_att = pd.concat(
        [
            pd.read_csv(
                os.path.join(tim_subj_dir, item),
                names=[item.split("_")[3].split(".")[0]],
            )
            for item in tf_list
        ],
        axis=1,
    )
    column_list = list(df_att)

    # new attribute (att) column, fill with column name,
    #   remove rows/vols w/multiple behaviors, fill
    #   NaN with "base"
    for beh in tim_beh_list:
        df_att.loc[df_att[beh] == 1, "att"] = beh

    for i, j in enumerate(df_att["att"]):
        if df_att[column_list].iloc[i].sum() > 1:
            df_att.at[i, "att"] = "base"
    df_att = df_att.replace(np.nan, "base", regex=True)

    # add category (cat) column for afni, 0 = base
    if "base" not in tim_beh_list:
        tim_beh_list.insert(0, "base")
    for cat, beh in enumerate(tim_beh_list):
        df_att.loc[df_att["att"] == beh, "cat"] = int(cat)
    df_att["cat"] = df_att["cat"].astype(int)

    # write
    h_out = os.path.join(tim_subj_dir, f"{tim_dcn_str}_categories.txt")
    np.savetxt(h_out, df_att["cat"].values, fmt="%s", delimiter=" ")

    h_df = os.path.join(tim_subj_dir, f"{tim_dcn_str}_cat_df.txt")
    np.savetxt(h_df, df_att.values, fmt="%s", delimiter=" ")


# %%
# For Testing
# work_dir = "/scratch/madlab/nate_vCAT/derivatives"
subj_dir = "/scratch/madlab/nate_vCAT/derivatives/sub-005/ses-S1"
# task_dict = {
#     "loc": ["face", "scene", "num"],
#     "Study": {"BE": [["Bfe", "Bse"], ["face", "scene"]]},
# }
task_dict = {"loc": ["face", "scene", "num"], "Study": {"BE": ["Bfe", "Bse"]}}
# phase = "loc"
# hdr_dict = func_lenRun(subj_dir, phase)

# %%
# make category files
for phase in task_dict:
    if type(task_dict[phase]) == list:
        hdr_dict = func_lenRun(
            subj_dir,
            phase,
        )
        func_timing(
            task_dict[phase],
            hdr_dict,
            subj_dir,
            phase,
            f"{phase}_decon",
        )
    elif type(task_dict[phase]) == dict:
        for decon in task_dict[phase]:
            hdr_dict = func_lenRun(
                subj_dir,
                phase,
            )
            func_timing(
                task_dict[phase][decon],
                hdr_dict,
                subj_dir,
                phase,
                f"{phase}_{decon}",
            )

# %%
