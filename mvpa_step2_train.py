"""
Notes

PyMVPA is written in Python 2.7
"""

# %%
import mvpa2
from mvpa2.suite import *
import os

# %%
"""
Example
"""
# ds = load_example_fmri_dataset()
# ds = ds[ds.chunks < 5]
# print ds.summary()

# # %%
# pl.figure(figsize=(14, 14))  # larger figure
# hist(ds, xgroup_attr='chunks', ygroup_attr='targets', noticks=None,
#      bins=20, normed=True)

# # next only works with floating point data
# ds.samples = ds.samples.astype('float')

# # look at sample similarity
# # Note, the decreasing similarity with increasing temporal distance
# # of the samples
# pl.figure(figsize=(14, 6))
# pl.subplot(121)
# plot_samples_distance(ds, sortbyattr='chunks')
# pl.title('Sample distances (sorted by chunks)')

# # %%
# # similar distance plot, but now samples sorted by their
# # respective targets, i.e. samples with same targets are plotted
# # in adjacent columns/rows.
# # Note, that the first and largest group corresponds to the
# # 'rest' condition in the dataset
# pl.subplot(122)
# plot_samples_distance(ds, sortbyattr='targets')
# pl.title('Sample distances (sorted by targets)')

# # z-score features individually per chunk
# print 'Detrending data'
# poly_detrend(ds, polyord=2, chunks_attr='chunks')
# print 'Z-Scoring data'
# zscore(ds)

# pl.figure(figsize=(14, 6))
# pl.subplot(121)
# plot_samples_distance(ds, sortbyattr='chunks')
# pl.title('Distances: z-scored, detrended (sorted by chunks)')
# pl.subplot(122)
# plot_samples_distance(ds, sortbyattr='targets')
# pl.title('Distances: z-scored, detrended (sorted by targets)')
# %%


"""
Data
"""
# set up
data_path = "/scratch/madlab/nate_vCAT/derivatives/mvpa"
mask_fname = os.path.join(
    data_path, "sub005/masks/orig/GM_int_mask.nii.gz")

dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
# dhandle = OpenFMRIDataset(data_path)
dhandle.get_subj_ids()
dhandle.get_task_descriptions()

# load data
task = 1
model = 1
subj = 5
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
real_ds = vstack(run_datasets, a=0)
print real_ds.summary()

# remove base
fds = fds[fds.sa.targets != 'base']
print fds.shape
