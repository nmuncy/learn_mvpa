"""
Notes

Written for Python 2.7 <-------------!!!
"""

# %%
from mvpa2.suite import *
import os

# %%
"""
Import data

Just one subject, for testing
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
fds_all = vstack(run_datasets, a=0)
print fds_all.summary()

# behavior functional dataset - remove base volumes
fds_beh = fds_all[fds_all.sa.targets != 'base']
# print fds_beh.shape
print fds_beh.summary()

# %%
"""
Plot
"""
# hist of feature values
pl.figure(figsize=(14, 14))  # larger figure
hist(fds_beh,
     xgroup_attr='chunks',
     ygroup_attr='targets',
     noticks=None,
     bins=20,
     normed=True)

# distance measure, sort by chunks or targets
fds_beh.samples = fds_beh.samples.astype('float')

pl.figure(figsize=(14, 6))
pl.subplot(121)
plot_samples_distance(fds_beh, sortbyattr='chunks')
pl.title('Sample distances (sorted by chunks)')

pl.figure(figsize=(14, 6))
pl.subplot(122)
plot_samples_distance(fds_beh, sortbyattr='targets')
pl.title('Sample distances (sorted by targets)')

# %%
"""
Train
"""
# subset dataset for face vs scene
fds_train = fds_beh[np.array(
    [l in ['face', 'scene'] for l in fds_beh.sa.targets], dtype='bool')]

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

# %%
"""
Sensitivity

This does not seem to agree with above results
"""
Clf = LinearCSVMC
svdmapper = SVDMapper()
get_SVD_sliced = lambda x: ChainMapper([svdmapper, StaticFeatureSelection(x)])

clfs = [('All orig.\nfeatures (%i)' % fds_train.nfeatures, Clf()),
        ('All Comps\n(%i)' % (fds_train.nsamples \
                 - (fds_train.nsamples / len(fds_train.UC)),),
                        MappedClassifier(Clf(), svdmapper)),
        ('First 30\nComp.', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(0, 30)))),
        ('Comp.\n31-60', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(31, 60)))),
        ('Comp.\n61-90', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(61, 90)))),
        ('Comp.\n91-120', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(91, 120)))),
        ('Comp.\n121-150', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(121, 150)))),
        ('Comp.\n151-180', MappedClassifier(Clf(),
                        get_SVD_sliced(slice(151, 180))))]

# run and visualize in barplot
results = []
labels = []

for desc, clf in clfs:
    print desc.replace('\n', ' ')
    # cv = CrossValidation(clf, NFoldPartitioner())
    cv = CrossValidation(clf,
                         NFoldPartitioner(),
                         errorfx=lambda p, t: np.mean(p == t))
    res = cv(fds_train)
    results.append(res.samples[:, 0])
    labels.append(desc)

plot_bars(results,
          labels=labels,
          title='Linear C-SVM classification (face vs. scene)',
          ylabel='Mean classification error (N-1 cross-validation, 12-fold)',
          distance=0.5)

# %%
"""
Load all
"""

# set up
data_path = "/scratch/madlab/nate_vCAT/derivatives/mvpa"
hdf5_path = os.path.join(data_path, "hdf5")
dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
dhandle.get_subj_ids()
dhandle.get_task_descriptions()

task = 1
model = 1
data_dict = {}

if not os.path.exists(os.path.join(hdf5_path, "data_all.hdf5.gz")):

    # load data
    for subj in dhandle.get_task_bold_run_ids(task):
        print subj

        run_datasets = []
        for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
            print run_id

            # pad, mask
            if subj < 10:
                subj_num = "sub00" + str(subj)
            else:
                subj_num = "sub0" + str(subj)

            mask_fname = os.path.join(data_path, subj_num,
                                      "masks/orig/GM_int_mask.nii.gz")
            run_events = dhandle.get_bold_run_model(model, subj, run_id)
            run_ds = dhandle.get_bold_run_dataset(subj,
                                                  task,
                                                  run_id,
                                                  chunks=run_id - 1,
                                                  mask=mask_fname)
            run_ds.sa['targets'] = events2sample_attr(run_events,
                                                      run_ds.sa.time_coords,
                                                      noinfolabel='base')
            run_datasets.append(run_ds)

        data_dict[subj] = run_datasets

    # save all data as one file
    mvpa2.base.hdf5.h5save(os.path.join(hdf5_path, "data_all.hdf5.gz"),
                           data_dict,
                           mode='w',
                           mkdir=True,
                           compression="gzip")

# %%
"""
Hyperalignment
"""

# example - get data
filepath = os.path.join(pymvpa_datadbroot,
                        'hyperalignment_tutorial_data.hdf5.gz')
ds_all = h5load(filepath)

# real - get data
fds_all = h5load(os.path.join(hdf5_path, "data_all.hdf5.gz"))

# example - zscore, inject the subject ID into all datasets
_ = [zscore(ds) for ds in ds_all]
# inject the subject ID into all datasets
for i, sd in enumerate(ds_all):
    sd.sa['subject'] = np.repeat(i, len(sd))

# %%
