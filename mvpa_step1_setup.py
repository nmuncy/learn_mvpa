"""
Notes
"""

# %%
import os
import sys
import json
import fnmatch
import subprocess
import pandas as pd
import numpy as np
from gp_step0_dcm2nii import func_sbatch


def func_lenRun(lrn_subj_dir, lrn_phase):

    """
    Returns a dictionary that contains TR length,
        number of volumes in the first run, number
        of runs in phase, and the length (sec) of a run.
    """

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


def func_timing(
    tim_beh_list,
    tim_hdr_dict,
    tim_subj_dir,
    tim_phase,
    tim_dcn_str,
):
    """
    Mines timing files (tf_phase_foo.txt) to construct
        (a) category file, used for training, and
        (b) matrix, used for verification.

    Files have one row per volume of phase, and category ints
        are 1-indexed from task_dict behaviors (0 = baseline)
    """

    # convert timing files to 1D in volume time
    for beh in tim_beh_list:

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
    #   update - censor 0 via 9999
    if "base" not in tim_beh_list:
        tim_beh_list.insert(0, "base")
    for cat, beh in enumerate(tim_beh_list):
        df_att.loc[df_att["att"] == beh, "cat"] = cat
    # df_att["cat"] = df_att["cat"].replace([0], 9999)
    df_att["cat"] = df_att["cat"].astype(int)

    # write categories, and matrix (for checking)
    h_out = os.path.join(tim_subj_dir, f"3dSVM_{tim_dcn_str}_categories.txt")
    np.savetxt(h_out, df_att["cat"].values, fmt="%s", delimiter=" ")
    h_df = os.path.join(tim_subj_dir, f"3dSVM_{tim_dcn_str}_matrix.txt")
    np.savetxt(h_df, df_att.values, fmt="%s", delimiter=" ")


def func_detrend(dtr_subj_dir, dtr_dcn_str, dtr_beh_list, dtr_hdr_dict):

    """
    Extracts non-baseline sub-bricks from decon_cbucket_REML+tlrc to
        generate 3dSVM_phase_foo+tlrc
    """

    # get relevant column numbers
    #   read design matrix, find column labels,
    #   then look for column label in task_dict
    #   to get a brick_list of only effects of interest
    brick_list = []
    with open(os.path.join(dtr_subj_dir, f"X.{dtr_dcn_str}.xmat.1D")) as f:
        h_file = f.readlines()
        for line in h_file:
            if line.__contains__("ColumnLabels"):
                col_list = line.split('"')[1].replace(" ", "").split(";")
                for i, j in enumerate(col_list):
                    if j.split("#")[0] in dtr_beh_list:
                        brick_list.append(f"{str(i)}")

    # determine correct number of sub-bricks
    #   decon adds an extra for model
    h_cmd = f"module load afni-20.2.06 \n 3dinfo -nv {dtr_subj_dir}/{dtr_dcn_str}_cbucket_REML+tlrc"
    h_len = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    len_wrong = h_len.communicate()[0].decode("utf-8").strip()
    len_right = int(len_wrong) - 2

    # make detrended data
    #   only contains tent data, extra sub-bricks
    #   are left behind
    h_cmd = f"""
        cd {dtr_subj_dir}
        3dTcat -prefix tmp_{dtr_dcn_str}_cbucket -tr {dtr_hdr_dict["LenTR"]} "{dtr_dcn_str}_cbucket_REML+tlrc[0..{len_right}]"
        3dSynthesize -prefix 3dSVM_{dtr_dcn_str} -matrix X.{dtr_dcn_str}.xmat.1D \
            -cbucket tmp_{dtr_dcn_str}_cbucket+tlrc -select {" ".join(brick_list)} -cenfill nbhr
    """
    subj_num = dtr_subj_dir.split("-")[1].split("/")[0]
    func_sbatch(h_cmd, 1, 1, 1, f"{subj_num}all", dtr_subj_dir)


def main():

    """
    Make category files and detrended EPI data.

    Account for whether type(task_dict[phase]) == list (one decon per phase),
        or type(task_dict[phase]) == dict (multiple decons per phase).

    Output will be titled 3dSVM_foo
    """

    # # For Testing
    # work_dir = "/scratch/madlab/nate_vCAT/derivatives"
    # subj = "sub-005"
    # sess = "ses-S1"
    # task_dict = {"loc": ["face", "scene"], "Study": {"BE": ["Bfe", "Bse"]}}
    # phase = "loc"
    # hdr_dict = func_lenRun(subj_dir, phase)

    # Receive args
    subj = str(sys.argv[1])
    sess = str(sys.argv[2])
    work_dir = str(sys.argv[3])
    subj_dir = os.path.join(work_dir, subj, sess)

    # Get dict
    with open(os.path.join(subj_dir, "task_dict.json")) as json_file:
        task_dict = json.load(json_file)

    # Work
    for phase in task_dict:
        if type(task_dict[phase]) == list:

            # get header info
            hdr_dict = func_lenRun(
                subj_dir,
                phase,
            )

            # set out str
            dcn_str = f"{phase}_decon"

            # make cat files
            if not os.path.exists(
                os.path.join(subj_dir, f"3dSVM_{dcn_str}_categories.txt")
            ):
                func_timing(
                    task_dict[phase],
                    hdr_dict,
                    subj_dir,
                    phase,
                    dcn_str,
                )

            # make detrended data
            if not os.path.exists(os.path.join(subj_dir, f"3dSVM_{dcn_str}+tlrc.HEAD")):
                func_detrend(subj_dir, dcn_str, task_dict[phase], hdr_dict)

        elif type(task_dict[phase]) == dict:
            for decon in task_dict[phase]:

                hdr_dict = func_lenRun(
                    subj_dir,
                    phase,
                )
                dcn_str = f"{phase}_{decon}"

                if not os.path.exists(
                    os.path.join(subj_dir, f"3dSVM_{dcn_str}_categories.txt")
                ):
                    func_timing(
                        task_dict[phase][decon],
                        hdr_dict,
                        subj_dir,
                        phase,
                        dcn_str,
                    )

                if not os.path.exists(
                    os.path.join(subj_dir, f"3dSVM_{dcn_str}+tlrc.HEAD")
                ):
                    func_detrend(subj_dir, dcn_str, task_dict[phase][decon], hdr_dict)


if __name__ == "__main__":
    main()
