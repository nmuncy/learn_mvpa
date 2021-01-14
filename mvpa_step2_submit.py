import os
import subprocess
from datetime import datetime


def main():

    # set up
    sess_str = "ses-S1"
    train_str = "loc_decon"

    deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    code_dir = "/home/nmuncy/compute/learn_mvpa"
    group_dir = os.path.join(deriv_dir, "grpAnalysis")
    mask = os.path.join(group_dir, "Group_Int_Mask.nii.gz")

    # submit job
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, f'Slurm_out/MVPA2_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    # Set stdout/err file
    h_out = os.path.join(out_dir, "out_train.txt")
    h_err = os.path.join(out_dir, "err_train.txt")

    # submit command
    sbatch_job = f"""
        sbatch \
        -J "MVPA2trn" -t 6:00:00 --mem=1000 --ntasks-per-node=1 \
        -p IB_44C_512G  -o {h_out} -e {h_err} \
        --account iacc_madlab --qos pq_madlab \
        --wrap="~/miniconda3/bin/python {code_dir}/mvpa_step2_train.py \
            {sess_str} {train_str} {deriv_dir} {group_dir} {mask}"
    """
    sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_submit.communicate()[0].decode("utf-8")
    print(job_id)


if __name__ == "__main__":
    main()
