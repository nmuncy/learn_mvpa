"""
Notes

This is written in Python 3.8
"""

# %%
import os
import fnmatch
import subprocess
from datetime import datetime
import time
import json

# set up
deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
mvpa_dir = os.path.join(deriv_dir, "mvpa")
code_dir = "/home/nmuncy/compute/learn_mvpa"
task_dict = {"loc": ["face", "scene", "num"]}
sess = "ses-S1"
beh_dur = 1
decon_type = "2GAM"


# %%
def main():
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
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()
    if not os.path.exists(mvpa_dir):
        os.makedirs(mvpa_dir)

    # write task_dict.json to avoid awkward import
    with open(os.path.join(mvpa_dir, "task_dict.json"), "w") as outfile:
        json.dump(task_dict, outfile)

    # scan key
    ref_subj = os.path.join(deriv_dir, subj_list[0])
    h_cmd = f"""
        module load afni-20.2.06
        3dinfo -tr {os.path.join(ref_subj, sess)}/run-1_{list(task_dict.keys())[0]}_scale+tlrc
    """
    h_tr = subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE)
    h_len_tr = h_tr.communicate()[0].decode("utf-8").strip()
    len_tr = float(h_len_tr)
    with open(os.path.join(mvpa_dir, "scan_key.txt"), "w") as scan_key:
        scan_key.write(f"TR {len_tr}")

    # study key
    with open(os.path.join(mvpa_dir, "study_key.txt"), "w") as study_key:
        study_key.write(f"{list(task_dict.keys())[0]} task")

    # task, condition keys
    if not os.path.exists(os.path.join(mvpa_dir, "task_key.txt")):
        for count, task in enumerate(task_dict.keys()):

            # task key
            with open(os.path.join(mvpa_dir, "task_key.txt"), "a") as task_key:
                task_key.write(f"task00{count+1} {task}")

            # condition keys
            cond_list = task_dict[task]
            for cc, cond in enumerate(cond_list):
                key_dir = os.path.join(mvpa_dir, "models", f"model00{count+1}")
                if not os.path.exists(key_dir):
                    os.makedirs(key_dir)
                with open(os.path.join(key_dir, "condition_key.txt"), "a") as cond_key:
                    cond_key.write(f"task00{count+1} cond00{cc+1} {cond}\n")

    # %%
    """
    Step 2: Group Intersection Mask

    PyMVPA requires the same number of voxels for each
    participant. Rather than use a template brain mask,
    use a group gray matter intersection mask.
    """
    group_dir = os.path.join(deriv_dir, "grpAnalysis")
    if not os.path.exists(group_dir):
        os.makedirs(group_dir)

    if not os.path.exists(os.path.join(group_dir, "Group_Int_Mask.nii.gz")):
        mask_list = []
        for subj in subj_list:
            mask_file = os.path.join(deriv_dir, subj, "ses-S1/mask_epi_anat+tlrc.HEAD")
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

    # %%
    """
    Step 3: Submit job for each subject
    """
    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, f'Slurm_out/MVPA1_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    for subj in subj_list:

        # Set stdout/err file
        h_out = os.path.join(out_dir, f"out_{subj}.txt")
        h_err = os.path.join(out_dir, f"err_{subj}.txt")

        # submit command
        subj_dir = os.path.join(deriv_dir, subj, sess)
        sbatch_job = f"""
            sbatch \
            -J "PPI{subj.split("-")[1]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
            -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="~/miniconda3/bin/python {code_dir}/mvpa_step1_setup.py \
                {subj} {subj_dir} {decon_type} {len_tr} {beh_dur} {deriv_dir}"
        """

        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)
        time.sleep(1)


if __name__ == "__main__":
    main()