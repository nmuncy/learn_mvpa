"""
Notes

Written for Python 2.7 <-------------!!!
"""

# %%
import mvpa2
from mvpa2.suite import *
import os

# %%
"""
Import data

How to import entire dataset?
"""
# set up
data_path = "/scratch/madlab/nate_vCAT/derivatives/mvpa"
mask_fname = os.path.join(data_path, "sub005/masks/orig/GM_int_mask.nii.gz")

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
    run_ds = dhandle.get_bold_run_dataset(subj,
                                          task,
                                          run_id,
                                          chunks=run_id - 1,
                                          mask=mask_fname)

    # convert event to sample attribute, assign as target
    run_ds.sa['targets'] = events2sample_attr(run_events,
                                              run_ds.sa.time_coords,
                                              noinfolabel='base')

    # write list
    run_datasets.append(run_ds)

# merge datasets, a=0 means attributes for first
#   set should be used for all data
real_ds = vstack(run_datasets, a=0)
print real_ds.summary()

# remove base
fds = real_ds[real_ds.sa.targets != 'base']
print fds.shape

# %%
"""
Plot
"""
# hist of feature values
pl.figure(figsize=(14, 14))  # larger figure
hist(fds,
     xgroup_attr='chunks',
     ygroup_attr='targets',
     noticks=None,
     bins=20,
     normed=True)

# distance measure, sort by chunks or targets
fds.samples = fds.samples.astype('float')

pl.figure(figsize=(14, 6))
pl.subplot(121)
plot_samples_distance(fds, sortbyattr='chunks')
pl.title('Sample distances (sorted by chunks)')

pl.figure(figsize=(14, 6))
pl.subplot(122)
plot_samples_distance(fds, sortbyattr='targets')
pl.title('Sample distances (sorted by targets)')

# %%
"""
Train
"""
# subset dataset for face vs scene
fds_train = fds[np.array([l in ['face', 'scene'] for l in fds.sa.targets],
                         dtype='bool')]

# set classifier
# clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
clf = LinearCSVMC()

# train, print
cvte = CrossValidation(clf,
                       NFoldPartitioner(),
                       errorfx=lambda p, t: np.mean(p == t),
                       enable_ca=['stats'])
cv_results = cvte(fds_train)
print cvte.ca.stats.as_string(description=True)
print cvte.ca.stats.matrix