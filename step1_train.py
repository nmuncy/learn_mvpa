"""
Notes
"""

# %%
from mvpa2.tutorial_suite import *

# %%
# set up
data_path = "/Users/nmuncy/Projects/learn_mvpa/vCAT_data"
mask_fname = os.path.join(data_path, "sub006", "masks", "orig", "GMc.nii.gz")

dhandle = OpenFMRIDataset(data_path)
dhandle.get_subj_ids()
dhandle.get_task_descriptions()

# load data
task = 1
model = 1
subj = 6
run_datasets = []

for run_id in dhandle.get_task_bold_run_ids(task)[subj]:

    # design info
    run_events = dhandle.get_bold_run_model(model, subj, run_id)

    # fmri data, mask
    run_ds = dhandle.get_bold_run_dataset(
        subj, task, run_id, chunks=run_id - 1, mask=mask_fname)

    # convert event to sample attribute, assign as target
    run_ds.sa['targets'] = events2sample_attr(
        run_events, run_ds.sa.time_coords, noinfolabel='base')

    # write list
    run_datasets.append(run_ds)

# merge datasets, a=0 means attributes for first
#   set should be used for all data
fds = vstack(run_datasets, a=0)
print fds.summary()

# remove base
fds = fds[fds.sa.targets != 'base']
print fds.shape

# %%
# Train
run_averager = mean_group_sample(['targets', 'chunks'])
fds = fds.get_mapped(run_averager)
fds.shape

# different cross-validation approach, get accuracy
clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
# clf = LinearCSVMC()
cvte = CrossValidation(clf, NFoldPartitioner(
), errorfx=lambda p, t: np.mean(p == t), enable_ca=['stats'])
cv_results = cvte(fds)
print cvte.ca.stats.as_string(description=True)
print cvte.ca.stats.matrix

# %%
