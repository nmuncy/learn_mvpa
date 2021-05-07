"""
Notes:

Wrapper script for step1_preproc.py.

Usage - update "set up" section.

phase_list = list of phases gathered within a single session.
    For example, if a study and then a test phase were both scanned
    during the same session, then phase_list = ["study", "test"]
"""
# %%
import os
from datetime import datetime
import subprocess
import time
import fnmatch

# set up
code_dir = "/home/nmuncy/compute/learn_mvpa"
parent_dir = "/scratch/madlab/nate_vCAT"
sess_dict = {"ses-S1": ["loc", "Study"]}


# %%
def main():

    # set up stdout/err capture
    current_time = datetime.now()
    out_dir = os.path.join(
        parent_dir,
        f'derivatives/Slurm_out/TS1_{current_time.strftime("%H%M_%d-%m-%y")}',
    )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # submit job for each subj/sess/phase
    subj_list = [
        x
        for x in os.listdir(os.path.join(parent_dir, "dset"))
        if fnmatch.fnmatch(x, "sub-*")
    ]
    subj_list.sort()

    for subj in subj_list:
        for sess in sess_dict:
            if not os.path.exists(
                os.path.join(
                    parent_dir,
                    "derivatives",
                    subj,
                    sess,
                    f"run-1_{sess_dict[sess][0]}_scale+tlrc.HEAD",
                )
            ):

                h_out = os.path.join(out_dir, f"out_{subj}_{sess}.txt")
                h_err = os.path.join(out_dir, f"err_{subj}_{sess}.txt")

                sbatch_job = f"""
                    sbatch \
                        -J "GP1{subj.split("-")[1]}" -t 10:00:00 --mem=4000 --ntasks-per-node=1 \
                        -p IB_44C_512G  -o {h_out} -e {h_err} \
                        --account iacc_madlab --qos pq_madlab \
                        --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
                        python {code_dir}/gp_step1_preproc.py {subj} {sess} \
                        {parent_dir} {' '.join(sess_dict[sess])}"
                """
                sbatch_submit = subprocess.Popen(
                    sbatch_job, shell=True, stdout=subprocess.PIPE
                )
                job_id = sbatch_submit.communicate()[0]
                print(job_id.decode("utf-8"))
                time.sleep(1)


if __name__ == "__main__":
    main()
