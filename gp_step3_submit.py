"""
Notes:

Wrapper script for step3_decon.py.

Update paths in "set up" section.

decon_type can be dmBLOCK, GAM, 2GAM, or TENT
"""

# %%
import os
from datetime import datetime
import fnmatch
import subprocess
import time
import json


# set up
code_dir = "/home/nmuncy/compute/learn_mvpa"
work_dir = "/scratch/madlab/nate_vCAT"
sess_list = ["ses-S1"]
# phase_list = ["loc", "Study"]
decon_type = "TENT"
decon_dict = {
    "loc": [
        "tf_loc_face.txt",
        "tf_loc_num.txt",
        "tf_loc_scene.txt",
    ],
    "Study": {
        "BP": [
            "tf_Study_Bfpc.txt",
            "tf_Study_Bfpi.txt",
            "tf_Study_Bspc.txt",
            "tf_Study_Bspi.txt",
        ],
        "BE": [
            "tf_Study_Bfec.txt",
            "tf_Study_Bfei.txt",
            "tf_Study_Bsec.txt",
            "tf_Study_Bsei.txt",
        ],
        "CP": [
            "tf_Study_cfpc.txt",
            "tf_Study_cfpi.txt",
            "tf_Study_cspc.txt",
            "tf_Study_cspi.txt",
        ],
        "CE": [
            "tf_Study_cfec.txt",
            "tf_Study_cfei.txt",
            "tf_Study_csec.txt",
            "tf_Study_csei.txt",
        ],
    },
}


def main():

    # set up stdout/err capture
    deriv_dir = os.path.join(work_dir, "derivatives")
    current_time = datetime.now()
    out_dir = os.path.join(
        deriv_dir, f'Slurm_out/TS3_{current_time.strftime("%H%M_%d-%m-%y")}'
    )
    os.makedirs(out_dir)

    # submit job for each subj/sess/phase
    subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
    subj_list.sort()

    # for i in subj_list:
    i = subj_list[0]
    for j in sess_list:

        h_out = os.path.join(out_dir, f"out_{i}_{j}.txt")
        h_err = os.path.join(out_dir, f"err_{i}_{j}.txt")

        # write decon_dict to json in subj dir
        with open(os.path.join(deriv_dir, i, j, "decon_dict.json"), "w") as outfile:
            json.dump(decon_dict, outfile)

        sbatch_job = f"""
            sbatch \
            -J "GP3{i.split("-")[1]}" -t 10:00:00 --mem=4000 --ntasks-per-node=1 \
            -p centos7_IB_44C_512G  -o {h_out} -e {h_err} \
            --account iacc_madlab --qos pq_madlab \
            --wrap="module load python-3.7.0-gcc-8.2.0-joh2xyk \n \
            python {code_dir}/gp_step3_decon.py {i} {j} {decon_type} {deriv_dir}"
        """
        # {' '.join(decon_dict.keys())}

        sbatch_submit = subprocess.Popen(sbatch_job, shell=True, stdout=subprocess.PIPE)
        job_id = sbatch_submit.communicate()[0]
        print(job_id)
        time.sleep(5)


if __name__ == "__main__":
    main()

# %%
