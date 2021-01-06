"""
Notes
"""

# %%
import os
import time
import fnmatch
import subprocess
from datetime import datetime
from gp_step0_dcm2nii import func_sbatch


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
                msk_deriv_dir, subj, "ses-S1/mask_epi_anat+tlrc.HEAD"
            )
            if os.path.exists(mask_file):
                mask_list.append(mask_file.split(".")[0])

        # combine anat-epi intersection masks of all subjs
        if not os.path.exists(os.path.join(group_dir, "Group_epi_int.nii.gz")):
            h_cmd = f"""
                # module load afni-20.2.06
                3dMean -prefix {group_dir}/Group_epi_mean.nii.gz {" ".join(mask_list)}
                3dmask_tool -input {" ".join(mask_list)} -frac 1 -prefix {group_dir}/Group_epi_int.nii.gz
            """
            # subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE).wait()
            func_sbatch(h_cmd, 1, 1, 1, "grpEpi", group_dir)

        # make GM mask, GM intersection mask
        if not os.path.exists(os.path.join(group_dir, "Group_Int_Mask.nii.gz")):
            atropos_dir = "/home/data/madlab/atlases/vold2_mni/priors_ACT"
            h_cmd = f"""
                module load c3d-1.0.0-gcc-8.2.0
                # module load afni-20.2.06

                cd {group_dir}
                c3d {atropos_dir}/Prior2.nii.gz {atropos_dir}/Prior4.nii.gz -add -o tmp_Prior_GM.nii.gz
                3dresample -master {mask_list[0]} -rmode NN -input tmp_Prior_GM.nii.gz -prefix tmp_Template_GM_mask.nii.gz

                c3d tmp_Template_GM_mask.nii.gz Group_epi_int.nii.gz -multiply -o tmp_Intersection_GM_prob_mask.nii.gz
                c3d tmp_Intersection_GM_prob_mask.nii.gz -thresh 0.1 1 1 0 -o Group_Int_Mask.nii.gz
                rm tmp_*
            """
            func_sbatch(h_cmd, 1, 1, 1, "grpInt", group_dir)
            # subprocess.Popen(h_cmd, shell=True, stdout=subprocess.PIPE).wait()


def func_submit(
    sub_code_dir, sub_deriv_dir, sub_subj_list, sub_sess, sub_model, sub_len_tr
):
    """
    Submit job for each subject
    """
    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        sub_deriv_dir, f'Slurm_out/SVM1_{current_time.strftime("%H%M_%d-%m-%y")}'
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
            -J "SVM1{subj.split("-")[1]}" -t 2:00:00 --mem=1000 --ntasks-per-node=1 \
            -p IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="~/miniconda3/bin/python {sub_code_dir}/3dsvm_step1_setup.py \
                {subj} {subj_dir} {sub_len_tr} {sub_deriv_dir} {sub_model}"
        """
        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)
        time.sleep(1)


# %%
def main():

    main_deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    main_task_dict = {
        "loc": ["face", "scene", "num"],
        "Study": {"BE": [["Bfe", "Bse"], ["face", "scene"]]},
    }

    # Make intersection mask
    subj_list = [x for x in os.listdir(main_deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()
    if not os.path.exists(
        os.path.join(main_deriv_dir, "grpAnalysis/Group_Int_Mask.nii.gz")
    ):
        func_mask(main_deriv_dir, subj_list)

    # submit job
