# %%
from mvpa2.tutorial_suite import *

# %%
"""
Searchlights
    We can classify, but where in the brain is
        the signal of interest?
"""

# load data
data_path = "/home/nate/Projects/pymvpa/tutorial_data_4/data"
mask_fname = os.path.join(data_path, "sub001", "masks", "orig", "vt.nii.gz")

dhandle = OpenFMRIDataset(data_path)
# dhandle.get_subj_ids()
# dhandle.get_task_descriptions()

task = 1
model = 1
subj = 1
run_datasets = []
for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
    run_events = dhandle.get_bold_run_model(model, subj, run_id)
    run_ds = dhandle.get_bold_run_dataset(
        subj, task, run_id, chunks=run_id - 1, mask=mask_fname)
    run_ds.sa['targets'] = events2sample_attr(
        run_events, run_ds.sa.time_coords, noinfolabel='rest')
    run_datasets.append(run_ds)

fds = vstack(run_datasets, a=0)
print fds.summary()

# preproc data
detrender = PolyDetrendMapper(polyord=1, chunks_attr='chunks')
detrended_fds = fds.get_mapped(detrender)
# print detrended_fds.a.mapper
zscore(detrended_fds, param_est=('targets', ['rest']))
fds = detrended_fds
# print fds.a.mapper
fds = fds[fds.sa.targets != 'rest']
print fds.shape

# set up for testing
rnames = {0: 'even', 1: 'odd'}
fds.sa['runtype'] = [rnames[c % 2] for c in fds.sa.chunks]
# averager = mean_group_sample(['targets', 'runtype'])
# type(averager)
# fds = fds.get_mapped(averager)
# fds.shape

# %%
# Stat test
aov = OneWayAnova()
f = aov(fds)
print f

# bound f-stat for fun
aov = OneWayAnova(postproc=FxMapper(
    'features', lambda x: x / x.max(), attrfx=None))
f = aov(fds)
print f.samples.max()

# %%
# set up cross-validation procedure, spheres
clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
cv = CrossValidation(clf, HalfPartitioner())
sl = sphere_searchlight(cv, radius=3, postproc=mean_sample())

# run on data
res = sl(fds)
print res

# do it with accuracy
cv = CrossValidation(clf, HalfPartitioner(
    attr='runtype'), errorfx=lambda p, t: np.mean(p == t))
sl = sphere_searchlight(cv, radius=3, postproc=mean_sample())
res = sl(fds)

sphere_acc = res.samples[0]
hist(sphere_acc, bins=np.linspace(0, 1, 18))
map2nifti(fds, 1.0 - sphere_acc).to_filename('sl.nii.gz')


# %%
"""
Repeat on all data (not just vt mask)
"""

run_datasets = []
for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
    run_events = dhandle.get_bold_run_model(model, subj, run_id)
    run_ds = dhandle.get_bold_run_dataset(
        subj, task, run_id, chunks=run_id - 1)
    run_ds.sa['targets'] = events2sample_attr(
        run_events, run_ds.sa.time_coords, noinfolabel='rest')
    run_datasets.append(run_ds)

fds = vstack(run_datasets, a=0)
print fds.summary()

# preproc data
detrender = PolyDetrendMapper(polyord=1, chunks_attr='chunks')
detrended_fds = fds.get_mapped(detrender)
# print detrended_fds.a.mapper
zscore(detrended_fds, param_est=('targets', ['rest']))
fds = detrended_fds
# print fds.a.mapper
fds = fds[fds.sa.targets != 'rest']
print fds.shape

# set up for testing
rnames = {0: 'even', 1: 'odd'}
fds.sa['runtype'] = [rnames[c % 2] for c in fds.sa.chunks]
averager = mean_group_sample(['targets', 'runtype'])
type(averager)
fds = fds.get_mapped(averager)
fds.shape

# set up cross-validation procedure, spheres
# clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
# cv = CrossValidation(clf, HalfPartitioner())
# sl = sphere_searchlight(cv, radius=3, postproc=mean_sample())

# run on data
res = sl(fds)
print res

sphere_errors = res.samples[0]
res_mean = np.mean(res)
res_std = np.std(res)
chance_level = 1.0 - (1.0 / len(fds.uniquetargets))
frac_lower = np.round(np.mean(sphere_errors < chance_level - 2 * res_std), 3)

hist(sphere_errors, bins=np.linspace(0, 1, 18))
map2nifti(fds, 1.0 - sphere_errors).to_filename('sl.nii.gz')
