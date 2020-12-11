"""
Notes

This is written in Python 3.8

Assumes gp_step1-3 have been run.

task_dict - for one deconvultion per phase, use list.
    for multiple decon per phase, use dict.

    Note - Test data needs to have same behaviors as Train, so
        task_dict[test][decon][0] = decon sub-brick labels
        task_dict[test][decon][1] = equivalent strings to study

Analyses are currently organized by models. Tasks are for separate
    deconvolutions.

    model001:
        train = loc (task001)
        test = Study BE (task002)
"""

import os
import fnmatch
import subprocess
from datetime import datetime
import time
import json


def func_keys(key_len_tr, key_study_str, key_mvpa_dir, key_task_dict, key_model):
    """
    Step 1: Keys

    1) Various keys are made so pymvpa can orient itself:

        scan_key.txt - contains TR length
        study_key.txt - has name of task (same string as in timing files)
        task_key.txt - identifies task number
        condition_key.txt identifies conditions of each task

    2) Pymvpa expects above keys in certain locations, and uses
        pseudo BIDS formatting for directory hierarchy, file names,
        and info in files.
    """

    # scan key - assumes all runs have same TR
    with open(os.path.join(key_mvpa_dir, "scan_key.txt"), "w") as scan_key:
        scan_key.write(f"TR {key_len_tr}")

    # study key
    with open(os.path.join(key_mvpa_dir, "study_key.txt"), "w") as study_key:
        study_key.write(key_study_str)

    # task key - account for single/multiple decons per phase
    if not os.path.exists(os.path.join(key_mvpa_dir, "task_key.txt")):
        for count, phase in enumerate(key_task_dict):
            with open(os.path.join(key_mvpa_dir, "task_key.txt"), "a") as task_key:
                if type(key_task_dict[phase]) == list:
                    task_key.write(f"task00{count+1} {phase}\n")
                elif type(key_task_dict[phase]) == dict:
                    for decon in key_task_dict[phase]:
                        task_key.write(f"task00{count+1} {decon}\n")

    # condition keys
    key_model_dir = os.path.join(key_mvpa_dir, f"models/{key_model}")
    if not os.path.exists(key_model_dir):
        os.makedirs(key_model_dir)

    cond_out = os.path.join(key_model_dir, "condition_key.txt")
    if not os.path.exists(cond_out):
        decon_count = 1
        for phase in key_task_dict:
            if type(key_task_dict[phase]) == list:
                for cc, cond in enumerate(key_task_dict[phase]):
                    with open(cond_out, "a") as cond_key:
                        cond_key.write(f"task00{decon_count} cond00{cc+1} {cond}\n")
                decon_count += 1
            elif type(key_task_dict[phase]) == dict:
                for decon in key_task_dict[phase]:
                    for cc, cond in enumerate(key_task_dict[phase][decon][1]):
                        with open(cond_out, "a") as cond_key:
                            cond_key.write(f"task00{decon_count} cond00{cc+1} {cond}\n")
                    decon_count += 1


def func_mask(msk_deriv_dir, msk_subj_list):
    """
    Step 2: Group Intersection Mask

    PyMVPA requires the same number of voxels for each
    participant. Rather than use a template brain mask,
    use a group gray matter intersection mask.

    Will be located in derivatives/grpAnalysis.
    """
    group_dir = os.path.join(msk_deriv_dir, "grpAnalysis")
    if not os.path.exists(group_dir):
        os.makedirs(group_dir)

    if not os.path.exists(os.path.join(group_dir, "Group_Int_Mask.nii.gz")):
        mask_list = []
        for subj in msk_subj_list:
            mask_file = os.path.join(
                msk_deriv_dir, subj, "ses-S1/mask_epi_anat+tlrc.HEAD"
            )
            if os.path.exists(mask_file):
                mask_list.append(mask_file.split(".")[0])

        # combine anat-epi intersection masks of all subjs
        if not os.path.exists(os.path.join(group_dir, "Group_epi_int.nii.gz")):
            h_cmd = f"""
                module load afni-20.2.06
                3dMean -prefix {group_dir}/Group_epi_mean.nii.gz {" ".join(mask_list)}
                3dmask_tool -input {" ".join(mask_list)} -frac 1 -prefix {group_dir}/Group_epi_int.nii.gz
            """
            subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE).wait()

        # make GM mask, GM intersection mask
        atropos_dir = "/home/data/madlab/atlases/vold2_mni/priors_ACT"
        h_cmd = f"""
            module load c3d/1.0.0
            module load afni-20.2.06

            cd {group_dir}
            c3d {atropos_dir}/Prior2.nii.gz {atropos_dir}/Prior4.nii.gz -add -o tmp_Prior_GM.nii.gz
            3dresample -master {mask_list[0]} -rmode NN -input tmp_Prior_GM.nii.gz -prefix tmp_Template_GM_mask.nii.gz

            c3d tmp_Template_GM_mask.nii.gz Group_epi_int.nii.gz -multiply -o tmp_Intersection_GM_prob_mask.nii.gz
            c3d tmp_Intersection_GM_prob_mask.nii.gz -thresh 0.1 1 1 0 -o Group_Int_Mask.nii.gz
            rm tmp_*
        """
        subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE).wait()


def func_submit(
    sub_code_dir, sub_deriv_dir, sub_subj_list, sub_sess, sub_model, sub_len_tr
):
    """
    Step 3: Submit job for each subject
    """
    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        sub_deriv_dir, f'Slurm_out/MVPA1_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    for subj in sub_subj_list:
        # subj = sub_subj_list[0]

        # Set stdout/err file
        h_out = os.path.join(out_dir, f"out_{subj}.txt")
        h_err = os.path.join(out_dir, f"err_{subj}.txt")

        # submit command
        subj_dir = os.path.join(sub_deriv_dir, subj, sub_sess)
        sbatch_job = f"""
            sbatch \
            -J "MVPA1{subj.split("-")[1]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
            -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="~/miniconda3/bin/python {sub_code_dir}/mvpa_step1_setup.py \
                {subj} {subj_dir} {sub_len_tr} {sub_deriv_dir} {sub_model}"
        """

        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)
        time.sleep(1)


def main():

    """
    Set up
    """
    main_deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    main_code_dir = "/home/nmuncy/compute/learn_mvpa"
    main_task_dict = {
        "loc": ["face", "scene", "num"],
        "Study": {"BE": [["Bfe", "Bse"], ["face", "scene"]]},
    }
    main_sess = "ses-S1"
    main_model = "model001"
    main_study = "vCAT Experiment"

    # get list of subjects
    main_subj_list = [
        x for x in os.listdir(main_deriv_dir) if fnmatch.fnmatch(x, "sub-*")
    ]
    main_subj_list.sort()

    # determine TR
    main_ref_file = os.path.join(
        main_deriv_dir,
        main_subj_list[0],
        main_sess,
        f"run-1_{list(main_task_dict.keys())[0]}_scale",
    )
    h_cmd = f"""
        module load afni-20.2.06
        3dinfo -tr {main_ref_file}
    """
    h_sub = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_len_tr = h_sub.communicate()[0].decode("utf-8").strip()
    main_len_tr = float(h_len_tr)

    # set, make output dirs, write task_dict
    main_mvpa_dir = os.path.join(main_deriv_dir, "mvpa")
    if not os.path.exists(main_mvpa_dir):
        os.makedirs(main_mvpa_dir)

    with open(os.path.join(main_mvpa_dir, "main_task_dict.json"), "w") as outfile:
        json.dump(main_task_dict, outfile)

    """
    Do work
    """
    # make keys
    func_keys(main_len_tr, main_study, main_mvpa_dir, main_task_dict, main_model)

    # make group mask
    func_mask(main_deriv_dir, main_subj_list)

    # submit sbatch jobs
    func_submit(
        main_code_dir,
        main_deriv_dir,
        main_subj_list,
        main_sess,
        main_model,
        main_len_tr,
    )


if __name__ == "__main__":
    main()
