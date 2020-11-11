"""
Notes

Written for Python 2.7 <-------------!!!
"""

# %%
from mvpa2.suite import *
import os
import gc

# %%
"""
Test run - do for one subj
"""

# # ## Import data
# # set up
# data_path = "/scratch/madlab/nate_vCAT/derivatives/mvpa"
# mask_fname = os.path.join(data_path, "sub005/masks/orig/GM_int_mask.nii.gz")

# dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
# # dhandle = OpenFMRIDataset(data_path)
# dhandle.get_subj_ids()
# dhandle.get_task_descriptions()

# # load data
# task = 1
# model = 1
# subj = 5
# run_datasets = []

# for run_id in dhandle.get_task_bold_run_ids(task)[subj]:

#     # design info
#     run_events = dhandle.get_bold_run_model(model, subj, run_id)

#     # fmri data, mask
#     run_ds = dhandle.get_bold_run_dataset(subj,
#                                           task,
#                                           run_id,
#                                           chunks=run_id - 1,
#                                           mask=mask_fname)

#     # convert event to sample attribute, assign as target
#     run_ds.sa['targets'] = events2sample_attr(run_events,
#                                               run_ds.sa.time_coords,
#                                               noinfolabel='base')

#     # write list
#     run_datasets.append(run_ds)

# # merge datasets, a=0 means attributes for first
# #   set should be used for all data
# fds_all = vstack(run_datasets, a=0)
# print fds_all.summary()

# # behavior functional dataset - remove base volumes
# fds_beh = fds_all[fds_all.sa.targets != 'base']
# # print fds_beh.shape
# print fds_beh.summary()

# # %%
# # ##Plot
# # hist of feature values
# pl.figure(figsize=(14, 14))  # larger figure
# hist(fds_beh,
#      xgroup_attr='chunks',
#      ygroup_attr='targets',
#      noticks=None,
#      bins=20,
#      normed=True)

# # distance measure, sort by chunks or targets
# fds_beh.samples = fds_beh.samples.astype('float')

# pl.figure(figsize=(14, 6))
# pl.subplot(121)
# plot_samples_distance(fds_beh, sortbyattr='chunks')
# pl.title('Sample distances (sorted by chunks)')

# pl.figure(figsize=(14, 6))
# pl.subplot(122)
# plot_samples_distance(fds_beh, sortbyattr='targets')
# pl.title('Sample distances (sorted by targets)')

# # %%
# # ##Train
# # subset dataset for face vs scene
# fds_train = fds_beh[np.array(
#     [l in ['face', 'scene'] for l in fds_beh.sa.targets], dtype='bool')]

# # set classifier
# # clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
# clf = LinearCSVMC()

# # train, print
# cvte = CrossValidation(clf,
#                        NFoldPartitioner(),
#                        errorfx=lambda p, t: np.mean(p == t),
#                        enable_ca=['stats'])
# cv_results = cvte(fds_train)
# print cvte.ca.stats.as_string(description=True)
# print cvte.ca.stats.matrix

# # %%
# # ##Sensitivity
# #   This does not seem to agree with above results

# Clf = LinearCSVMC
# svdmapper = SVDMapper()


# def get_SVD_sliced(x):
#     return ChainMapper([svdmapper, StaticFeatureSelection(x)])


# clfs = [
#     ('All orig.\nfeatures (%i)' % fds_train.nfeatures, Clf()),
#     ('All Comps\n(%i)' % (fds_train.nsamples -
#                           (fds_train.nsamples / len(fds_train.UC)), ),
#      MappedClassifier(Clf(), svdmapper)),
#     ('First 30\nComp.', MappedClassifier(Clf(), get_SVD_sliced(slice(0, 30)))),
#     ('Comp.\n31-60', MappedClassifier(Clf(), get_SVD_sliced(slice(31, 60)))),
#     ('Comp.\n61-90', MappedClassifier(Clf(), get_SVD_sliced(slice(61, 90)))),
#     ('Comp.\n91-120', MappedClassifier(Clf(), get_SVD_sliced(slice(91, 120)))),
#     ('Comp.\n121-150', MappedClassifier(Clf(), get_SVD_sliced(slice(121,
#                                                                     150)))),
#     ('Comp.\n151-180', MappedClassifier(Clf(), get_SVD_sliced(slice(151,
#                                                                     180))))
# ]

# # run and visualize in barplot
# results = []
# labels = []

# for desc, clf in clfs:
#     print desc.replace('\n', ' ')
#     # cv = CrossValidation(clf, NFoldPartitioner())
#     cv = CrossValidation(clf,
#                          NFoldPartitioner(),
#                          errorfx=lambda p, t: np.mean(p == t))
#     res = cv(fds_train)
#     results.append(res.samples[:, 0])
#     labels.append(desc)

# plot_bars(results,
#           labels=labels,
#           title='Linear C-SVM classification (face vs. scene)',
#           ylabel='Mean classification error (N-1 cross-validation, 12-fold)',
#           distance=0.5)

# %%
"""
Real Work

MVPA with all data

Step 0: Load data
"""

# set up
data_path = "/Users/nmuncy/Projects/learn_mvpa/mvpa"
tutorial_path = os.path.join(data_path, "datadb")
hdf5_path = os.path.join(data_path, "hdf5")

# %%
dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
dhandle.get_subj_ids()
dhandle.get_task_descriptions()

task = 1
model = 1
data_list = []

if not os.path.exists(os.path.join(hdf5_path, "data_all.hdf5.gz")):

    # load data
    for subj in dhandle.get_task_bold_run_ids(task):
        # subj = 5
        print subj

        run_datasets = []
        for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
            print run_id

            # pad
            if subj < 10:
                subj_num = "sub00" + str(subj)
            else:
                subj_num = "sub0" + str(subj)

            # determine mask
            mask_fname = os.path.join(data_path, subj_num,
                                      "masks/orig/Group_Int_Mask.nii.gz")

            # get model
            run_events = dhandle.get_bold_run_model(model, subj, run_id)

            # get data
            run_ds = dhandle.get_bold_run_dataset(subj,
                                                  task,
                                                  run_id,
                                                  chunks=run_id - 1,
                                                  mask=mask_fname)

            # get behavior events
            run_ds.sa['targets'] = events2sample_attr(run_events,
                                                      run_ds.sa.time_coords,
                                                      noinfolabel='base')

            run_datasets.append(run_ds)

        # check that number of voxels (length) is equal across subjs
        #   this is a result of the mask
        print run_ds.fa.values()

        # append list, collapse across runs
        data_list.append(vstack(run_datasets, a=0))

    # save all data as one file
    mvpa2.base.hdf5.h5save(os.path.join(hdf5_path, "data_all.hdf5.gz"),
                           data_list,
                           mode='w',
                           mkdir=True,
                           compression="gzip")

del run_ds
del data_list
del run_datasets

# %%
"""
Step 2: Between, Within Classifiers
"""

# # example - get data
# ds_all = h5load(
#     os.path.join(tutorial_path, "hyperalignment_tutorial_data.hdf5.gz"))

# real - get data
fds_all = h5load(os.path.join(hdf5_path, "data_all.hdf5.gz"))

# # example - zscore, inject the subject ID into all datasets
# _ = [zscore(ds) for ds in ds_all]
# for i, sd in enumerate(ds_all):
#     sd.sa['subject'] = np.repeat(i, len(sd))

# nsubjs = len(ds_all)
# ncats = len(ds_all[0].UT)
# nruns = len(ds_all[0].UC)

# real - match
for i, sd in enumerate(fds_all):
    sd.sa['subject'] = np.repeat(i, len(sd))
nsubjs_real = len(fds_all)
ncats_real = len(fds_all[0].UT)
nruns_real = len(fds_all[0].UC)

# %%
# Classifier, feature selection - 100 highest
clf = LinearCSVMC()
nf = 100

fselector = FixedNElementTailSelector(nf,
                                      tail='upper',
                                      mode='select',
                                      sort=False)

sbfs = SensitivityBasedFeatureSelection(OneWayAnova(),
                                        fselector,
                                        enable_ca=['sensitivities'])

fsclf = FeatureSelectionClassifier(clf, sbfs)

# Classify within, between
cvws = CrossValidation(fsclf,
                       NFoldPartitioner(attr='chunks'),
                       errorfx=mean_match_accuracy)

cvbs = CrossValidation(fsclf,
                       NFoldPartitioner(attr='subject'),
                       errorfx=mean_match_accuracy)

# # example - train within/between-subject classifier
# wsc_results = [cvws(sd) for sd in ds_all]
# wsc_results = vstack(wsc_results)
# ds_comb = vstack(ds_all)
# bsc_results = cvbs(ds_comb)

# real
wsc_real = [cvws(sd) for sd in fds_all]
wsc_real = vstack(wsc_real)
fds_comb = vstack(fds_all)
bsc_real = cvbs(fds_comb)

# # example - compare
# print "Example within-subject: %.2f +/-%.3f" % (
#     np.mean(wsc_results), np.std(wsc_results) / np.sqrt(nsubjs - 1))

# print "Example between-subject (anatomically aligned): %.2f +/-%.3f" % (np.mean(
#     bsc_results), np.std(np.mean(bsc_results, axis=1)) / np.sqrt(nsubjs - 1))

# real- compare
print "Average Classification Accuracies"

print "within-subject: %.2f +/-%.3f" % (
    np.mean(wsc_real), np.std(wsc_real) / np.sqrt(nsubjs_real - 1))

print "between-subject (anatomically aligned): %.2f +/-%.3f" % (np.mean(
    bsc_real), np.std(np.mean(bsc_real, axis=1)) / np.sqrt(nsubjs_real - 1))


# %%
"""
Step 3: Hyperalignment
"""
# classifier
cv = CrossValidation(clf, NFoldPartitioner(attr='subject'),
                     errorfx=mean_match_accuracy)

# # example
# bsc_hyper_results = []
# for test_run in range(nruns):
#     ds_train = [sd[sd.sa.chunks != test_run, :] for sd in ds_all]
#     ds_test = [sd[sd.sa.chunks == test_run, :] for sd in ds_all]

#     anova = OneWayAnova()
#     fscores = [anova(sd) for sd in ds_train]
#     featsels = [StaticFeatureSelection(fselector(fscore))
#                 for fscore in fscores]
#     ds_train_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_train)]

#     hyper = Hyperalignment()
#     hypmaps = hyper(ds_train_fs)

#     ds_test_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_test)]
#     ds_hyper = [h.forward(sd) for h, sd in zip(hypmaps, ds_test_fs)]

#     ds_hyper = vstack(ds_hyper)
#     zscore(ds_hyper, chunks_attr='subject')
#     res_cv = cv(ds_hyper)
#     bsc_hyper_results.append(res_cv)

# bsc_hyper_results = hstack(bsc_hyper_results)


# Real
bsc_hyper_real = []
for test_run in range(nruns_real):
    # test_run = 1
    ds_train = [sd[sd.sa.chunks != test_run, :] for sd in fds_all]
    ds_test = [sd[sd.sa.chunks == test_run, :] for sd in fds_all]

    anova = OneWayAnova()
    fscores = [anova(sd) for sd in ds_train]
    featsels = [StaticFeatureSelection(fselector(fscore))
                for fscore in fscores]
    ds_train_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_train)]

    hyper = Hyperalignment()
    hypmaps = hyper(ds_train_fs)

    ds_test_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_test)]
    ds_hyper = [h.forward(sd) for h, sd in zip(hypmaps, ds_test_fs)]

    ds_hyper = vstack(ds_hyper)
    zscore(ds_hyper, chunks_attr='subject')
    res_cv = cv(ds_hyper)
    bsc_hyper_real.append(res_cv)

bsc_hyper_real = hstack(bsc_hyper_real)


# # example - compare
# print "Example within-subject: %.2f +/-%.3f" % (
#     np.mean(wsc_results), np.std(wsc_results) / np.sqrt(nsubjs - 1))

# print "Example between-subject (anatomically aligned): %.2f +/-%.3f" % (np.mean(
#     bsc_results), np.std(np.mean(bsc_results, axis=1)) / np.sqrt(nsubjs - 1))

# print "Example between-subject (hyperaligned): %.2f +/-%.3f" % (np.mean(
#     bsc_hyper_results), np.std(np.mean(bsc_hyper_results, axis=1)) / np.sqrt(nsubjs - 1))


# real - compare
print "Average Classification Accuracies"

print "within-subject: %.2f +/-%.3f" % (
    np.mean(wsc_real), np.std(wsc_real) / np.sqrt(nsubjs_real - 1))

print "between-subject (anatomically aligned): %.2f +/-%.3f" % (np.mean(
    bsc_real), np.std(np.mean(bsc_real, axis=1)) / np.sqrt(nsubjs_real - 1))

print "between-subject (hyperaligned): %.2f +/-%.3f" \
    % (np.mean(bsc_hyper_real),
       np.std(np.mean(bsc_hyper_real, axis=1)) / np.sqrt(nsubjs_real - 1))


# # example - similarity
# anova = OneWayAnova()
# fscores = [anova(sd) for sd in ds_all]
# fscores = np.mean(np.asarray(vstack(fscores)), axis=0)
# ds_fs = [sd[:, fselector(fscores)] for sd in ds_all]
# hyper = Hyperalignment()
# mappers = hyper(ds_fs)
# ds_hyper = [m.forward(ds_) for m, ds_ in zip(mappers, ds_fs)]
# sm_orig = [np.corrcoef(
#     sd.get_mapped(
#         mean_group_sample(['targets'])).samples)
#     for sd in ds_fs]
# sm_orig_mean = np.mean(sm_orig, axis=0)

# sm_hyper_mean = np.mean(
#     [np.corrcoef(
#         sd.get_mapped(mean_group_sample(['targets'])).samples)
#      for sd in ds_hyper],
#     axis=0)
# ds_hyper = vstack(ds_hyper)
# sm_hyper = np.corrcoef(ds_hyper.get_mapped(mean_group_sample(['targets'])))
# ds_fs = vstack(ds_fs)
# sm_anat = np.corrcoef(ds_fs.get_mapped(mean_group_sample(['targets'])))

# intended_label_order = [2, 4, 1, 5, 3, 0, 6]
# labels = ds_all[0].UT
# labels = labels[intended_label_order]

# pl.figure(figsize=(12, 12))
# # plot all three similarity structures
# for i, sm_t in enumerate((
#     (sm_orig_mean, "Average within-subject\nsimilarity"),
#     (sm_anat, "Similarity of group average\ndata (anatomically aligned)"),
#     (sm_hyper_mean, "Average within-subject\nsimilarity (hyperaligned data)"),
#     (sm_hyper, "Similarity of group average\ndata (hyperaligned)"),
# )):
#     sm, title = sm_t
#     # reorder matrix columns to match label order
#     sm = sm[intended_label_order][:, intended_label_order]
#     pl.subplot(2, 2, i + 1)
#     pl.imshow(sm, vmin=-1.0, vmax=1.0, interpolation='nearest')
#     pl.colorbar(shrink=.4, ticks=[-1, 0, 1])
#     pl.title(title, size=12)
#     ylim = pl.ylim()
#     pl.xticks(range(ncats), labels, size='small', stretch='ultra-condensed',
#               rotation=45)
#     pl.yticks(range(ncats), labels, size='small', stretch='ultra-condensed',
#               rotation=45)
#     pl.ylim(ylim)


# real data - similarity
anova = OneWayAnova()
fscores = [anova(sd) for sd in fds_all]
fscores = np.mean(np.asarray(vstack(fscores)), axis=0)

ds_fs = [sd[:, fselector(fscores)] for sd in fds_all]
hyper = Hyperalignment()
mappers = hyper(ds_fs)
ds_hyper = [m.forward(ds_) for m, ds_ in zip(mappers, ds_fs)]

sm_orig = [np.corrcoef(
    sd.get_mapped(
        mean_group_sample(['targets'])).samples)
    for sd in ds_fs]
sm_orig_mean = np.mean(sm_orig, axis=0)

sm_hyper_mean = np.mean(
    [np.corrcoef(
        sd.get_mapped(mean_group_sample(['targets'])).samples)
     for sd in ds_hyper],
    axis=0)

ds_hyper = vstack(ds_hyper)
sm_hyper = np.corrcoef(ds_hyper.get_mapped(mean_group_sample(['targets'])))

ds_fs = vstack(ds_fs)
sm_anat = np.corrcoef(ds_fs.get_mapped(mean_group_sample(['targets'])))

intended_label_order = [0, 2, 1, 3]
labels = fds_all[0].UT
labels = labels[intended_label_order]

pl.figure(figsize=(12, 12))
# plot all three similarity structures
for i, sm_t in enumerate((
    (sm_orig_mean, "Average within-subject\nsimilarity"),
    (sm_anat, "Similarity of group average\ndata (anatomically aligned)"),
    (sm_hyper_mean, "Average within-subject\nsimilarity (hyperaligned data)"),
    (sm_hyper, "Similarity of group average\ndata (hyperaligned)"),
)):
    sm, title = sm_t
    # reorder matrix columns to match label order
    sm = sm[intended_label_order][:, intended_label_order]
    pl.subplot(2, 2, i + 1)
    pl.imshow(sm, vmin=-1.0, vmax=1.0, interpolation='nearest')
    pl.colorbar(shrink=.4, ticks=[-1, 0, 1])
    pl.title(title, size=12)
    ylim = pl.ylim()
    pl.xticks(range(ncats_real), labels, size='small', stretch='ultra-condensed',
              rotation=45)
    pl.yticks(range(ncats_real), labels, size='small', stretch='ultra-condensed',
              rotation=45)
    pl.ylim(ylim)
# %%
