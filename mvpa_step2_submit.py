"""
Notes

This is written in Python 2.7
"""

# %%
import os
import subprocess
from datetime import datetime
import time


def func_job(deriv_dir, code_dir):
    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir,
        'Slurm_out/MVPA2_{}'.format(current_time.strftime("%H%M_%d-%m-%y")))

    os.makedirs(out_dir)
    h_out = os.path.join(out_dir, "out.txt")
    h_err = os.path.join(out_dir, "err.txt")

    # submit command
    sbatch_job = """
        sbatch \
        -J "MVPA2" -t 10:00:00 --mem=4000 --ntasks-per-node=10 \
        -p centos7_IB_44C_512G  -o {} -e {} \
        --account iacc_madlab --qos pq_madlab \
        --wrap="/home/data/madlab/envs/mvpa_env/bin/python {}/mvpa_step2_train.py {}"
    """.format(h_out, h_err, code_dir, deriv_dir)

    sbatch_submit = subprocess.Popen(sbatch_job,
                                     shell=True,
                                     stdout=subprocess.PIPE)
    job_id = sbatch_submit.communicate()[0]
    print(job_id)


def main():
    h_deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
    h_code_dir = "/home/nmuncy/compute/learn_mvpa"
    func_job(h_deriv_dir, h_code_dir)


if __name__ == "__main__":
    main()
# %%
