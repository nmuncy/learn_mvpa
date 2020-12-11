"""
Notes

This is written in Python 2.7

Make sure mvpa_env is loaded
"""

# %%
import os
import subprocess
import json
from datetime import datetime


def main():
    """
    model001 = loc, BE data
        task001 = loc - training data
        task002 = BE - test data

    main_model_dict = {<model num>: {"Train": <task num>, "Test": <task num>}}
    """

    deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    code_dir = "/home/nmuncy/compute/learn_mvpa"
    model_dict = {1: {"Train": 1, "Test": 2}}

    # write decon_dict to json in subj dir
    with open(os.path.join(deriv_dir, "mvpa/model_dict.json"), "w") as outfile:
        json.dump(model_dict, outfile)

    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, "Slurm_out/MVPA2_{}".format(current_time.strftime("%H%M_%d-%m-%y"))
    )

    os.makedirs(out_dir)
    h_out = os.path.join(out_dir, "out.txt")
    h_err = os.path.join(out_dir, "err.txt")

    # submit command
    # --account iacc_madlab --qos pq_madlab \
    sbatch_job = """
        sbatch \
        -J "MVPA2" -t 10:00:00 --mem=4000 --ntasks-per-node=10 \
        -p centos7_IB_44C_512G  -o {} -e {} \
        --wrap="/home/data/madlab/envs/mvpa_env/bin/python {}/mvpa_step2_build.py {}"
    """.format(
        h_out, h_err, code_dir, deriv_dir
    )
    sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
    job_id = sbatch_submit.communicate()[0]
    print(job_id)


if __name__ == "__main__":
    main()
# %%
