"""
Notes

task_dict - training and test behaviors must be equal in size (e.g. 2 & 2)
    so, if "loc" is train data, and training on "face" and "scene",
    test phases (Study) must also have 2 behaviors.

    key: [list] => decon_string: [tf_strA, tf_strB]

TODO update to account for incorrect CE/P behaviors
"""

# %%
import os
import time
import json
import fnmatch
import subprocess
from datetime import datetime
from gp_step0_dcm2nii import func_sbatch


# %%
# set up
sess = "ses-S1"
deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
code_dir = "/home/nmuncy/compute/learn_mvpa"
task_dict = {
    "loc": ["face", "scene"],
    "Study": {
        "BE": ["Bfe", "Bse"],
        "CE": ["cfec", "csec"],
        "FP": ["Ffpc", "Fspc"],
    },
}
atropos_dir = "/home/data/madlab/atlases/vold2_mni/priors_ACT"


# %%
def func_mask(msk_deriv_dir, msk_subj_list):
    """
    Group Intersection Mask

    Rather than use a template brain mask,
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
                msk_deriv_dir, subj, f"{sess}/mask_epi_anat+tlrc.HEAD"
            )
            if os.path.exists(mask_file):
                mask_list.append(mask_file.split(".")[0])

        # combine anat-epi intersection masks of all subjs
        if not os.path.exists(os.path.join(group_dir, "Group_epi_int.nii.gz")):
            h_cmd = f"""
                3dMean -prefix {group_dir}/Group_epi_mean.nii.gz {" ".join(mask_list)}

                3dmask_tool \
                    -input {" ".join(mask_list)} \
                    -frac 1 \
                    -prefix {group_dir}/Group_epi_int.nii.gz
            """
            func_sbatch(h_cmd, 1, 1, 1, "grpEpi", group_dir)

        # make GM mask, GM intersection mask
        if not os.path.exists(os.path.join(group_dir, "Group_Int_Mask.nii.gz")):
            h_cmd = f"""
                module load c3d-1.0.0-gcc-8.2.0
                cd {group_dir}

                c3d \
                    {atropos_dir}/Prior2.nii.gz \
                    {atropos_dir}/Prior4.nii.gz \
                    -add \
                    -o tmp_Prior_GM.nii.gz

                3dresample \
                    -master {mask_list[0]} \
                    -rmode NN \
                    -input tmp_Prior_GM.nii.gz \
                    -prefix tmp_Template_GM_mask.nii.gz

                c3d \
                    tmp_Template_GM_mask.nii.gz \
                    Group_epi_int.nii.gz \
                    -multiply \
                    -o tmp_Intersection_GM_prob_mask.nii.gz

                c3d \
                    tmp_Intersection_GM_prob_mask.nii.gz \
                    -thresh 0.1 1 1 0 \
                    -o Group_Int_Mask.nii.gz

                rm tmp_*
            """
            func_sbatch(h_cmd, 1, 1, 1, "grpInt", group_dir)


# %%
def main():

    """ Step 1: Make Intersection Mask"""
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()
    if not os.path.exists(os.path.join(deriv_dir, "grpAnalysis/Group_Int_Mask.nii.gz")):
        func_mask(deriv_dir, subj_list)

    """ Step 2: Submit Jobs"""
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, f'Slurm_out/MVPA1_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    for subj in subj_list:
        # subj = subj_list[0]

        # write json
        subj_dir = os.path.join(deriv_dir, subj, sess)
        with open(os.path.join(subj_dir, "task_dict.json"), "w") as outfile:
            json.dump(task_dict, outfile)

        # Set stdout/err file
        h_out = os.path.join(out_dir, f"out_{subj}.txt")
        h_err = os.path.join(out_dir, f"err_{subj}.txt")

        # submit command
        sbatch_job = f"""
            sbatch \
            -J "MVPA1{subj.split("-")[1]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
            -p IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="~/miniconda3/bin/python {code_dir}/mvpa_step1_setup.py \
                {subj} {sess} {deriv_dir}"
        """
        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)
        time.sleep(1)


if __name__ == "__main__":
    main()
