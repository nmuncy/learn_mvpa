# %%
import os
import fnmatch
import pandas as pd
import numpy as np

# %%
deriv_dir = "/scratch/madlab/nate_vCAT/derivatives"
beh_dir = os.path.join(deriv_dir, "grpAnalysis/behAnalysis")
test_list = ["Study_BE", "Study_BP", "Study_CP", "Study_FP"]
sess_str = "ses-S1"
cat_dict = {"Face": 1, "Scene": 2}

if not os.path.exists(beh_dir):
    os.makedirs(beh_dir)

# %%

subj_list = [x for x in os.listdir(deriv_dir) if fnmatch.fnmatch(x, "sub-*")]
subj_list.sort

# for test in test_list:
test = test_list[0]

# start master df
master_df = pd.DataFrame()

# for subj in subj_list:
subj = subj_list[1]
subj_dir = os.path.join(deriv_dir, subj, sess_str)

# Combine true, predicted categories
df_subj = pd.read_table(os.path.join(subj_dir, f"3dSVM_{test}_cat_updated.txt"))
df_subj["Pred"] = pd.read_table(os.path.join(subj_dir, f"3dSVM_pred_{test}.1D"))
df_subj.columns = ["True", "Pred"]

# Split into individual cat files, for ROC
for cat in cat_dict:
    np.savetxt(
        os.path.join(subj_dir, f"3dSVM_table_{test}_{cat}.txt"),
        df_subj[df_subj["True"].astype(str).str.contains(f"{cat_dict[cat]}")],
        fmt="%s",
        delimiter=",",
        header="True,Pred",
        comments="",
    )


# %%
