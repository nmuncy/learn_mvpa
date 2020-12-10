"""
Notes

Written for Python 2.7 <-------------!!!

TODO:
    1) Train for specific categories in block set?
    2) Verify that I am using separate data for
        determining hyperalignment parameter and
        training classifier
    3) Train on loc phase, test on vCAT phase
"""

# %%
from mvpa2.suite import *
import os
import gc
import sys
import json


# %%
def func_one_subj(data_path):
    """
    Step 0: Test run for one subject

    """

    # ## Import data
    # set up
    mask_subj = os.path.join(data_path,
                             "sub005/masks/orig/Group_Int_Mask.nii.gz")

    dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
    dhandle.get_subj_ids()
    dhandle.get_task_descriptions()

    # load data
    task = 1
    model = 1
    subj = 5
    run_data = []

    for run_id in dhandle.get_task_bold_run_ids(task)[subj]:

        # design info
        run_events = dhandle.get_bold_run_model(model, subj, run_id)

        # fmri data, mask
        run_data = dhandle.get_bold_run_dataset(subj,
                                                task,
                                                run_id,
                                                chunks=run_id - 1,
                                                mask=mask_subj)

        # convert event to sample attribute, assign as target
        run_data.sa['targets'] = events2sample_attr(run_events,
                                                    run_data.sa.time_coords,
                                                    noinfolabel='base')

        # write list
        run_data.append(run_data)

    # merge datasets, a=0 means attributes for first
    #   set should be used for all data
    fds_train = vstack(run_data, a=0)
    print fds_train.summary()

    # behavior functional dataset - remove base volumes
    fds_beh = fds_train[fds_train.sa.targets != 'base']
    # print fds_beh.shape
    print fds_beh.summary()

    # Plot
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

    # Train
    #   subset dataset for face vs scene
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

    # Sensitivity
    #   This does not seem to agree with above results

    Clf = LinearCSVMC
    svdmapper = SVDMapper()

    def get_SVD_sliced(x):
        return ChainMapper([svdmapper, StaticFeatureSelection(x)])

    clfs = [('All orig.\nfeatures (%i)' % fds_train.nfeatures, Clf()),
            ('All Comps\n(%i)' % (fds_train.nsamples -
                                  (fds_train.nsamples / len(fds_train.UC)), ),
             MappedClassifier(Clf(), svdmapper)),
            ('First 30\nComp.',
             MappedClassifier(Clf(), get_SVD_sliced(slice(0, 30)))),
            ('Comp.\n31-60',
             MappedClassifier(Clf(), get_SVD_sliced(slice(31, 60)))),
            ('Comp.\n61-90',
             MappedClassifier(Clf(), get_SVD_sliced(slice(61, 90)))),
            ('Comp.\n91-120',
             MappedClassifier(Clf(), get_SVD_sliced(slice(91, 120)))),
            ('Comp.\n121-150',
             MappedClassifier(Clf(), get_SVD_sliced(slice(121, 150)))),
            ('Comp.\n151-180',
             MappedClassifier(Clf(), get_SVD_sliced(slice(151, 180))))]

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

    plot_bars(
        results,
        labels=labels,
        title='Linear C-SVM classification (face vs. scene)',
        ylabel='Mean classification error (N-1 cross-validation, 12-fold)',
        distance=0.5)


# %%
def func_make_hdf5(mkhd_data_path, mkhd_hdf5_path, mkhd_model_dict):
    """
    Make training data (task001)
        Exclude base, num attributes (targets) from training

    Load all data in mvpa/sub*, save to
        mkhd_hdf5_path.

    Update - dhandle.get_bold_run_model() seems not to be able
        to handle different tasks, but stupidly loads all task
        data. Reverted to manual way of making data.
            Will this break the classifier?
    """
    # # For testing
    # main_par_path = "/Users/nmuncy/Projects/learn_mvpa"
    # mkhd_data_path = os.path.join(main_par_path, "mvpa")
    # mkhd_hdf5_path = os.path.join(mkhd_data_path, "hdf5")
    # mkhd_model_dict = {1: {"Train": 1, "Test": 2}}

    # make handler, load dataset info
    dhandle = mvpa2.datasets.sources.OpenFMRIDataset(mkhd_data_path)
    dhandle.get_subj_ids()
    dhandle.get_task_descriptions()

    # make train, test hdf5 files for each planned model
    for model in mkhd_model_dict:
        # model = 1
        print "model{}".format(model)

        for data_type in mkhd_model_dict[model]:
            # data_type = "Train"
            print data_type

            # set task, model, out_file string
            task = mkhd_model_dict[model][data_type]
            print "task{}".format(task)
            out_file = os.path.join(
                mkhd_hdf5_path, "model{}_data_{}.hdf5.gz".format(model, data_type))
            if not os.path.exists(out_file):

                # load data
                group_data = []
                for subj in dhandle.get_task_bold_run_ids(task):
                    # subj = 5
                    print "sub{}".format(subj)

                    subj_data = []
                    for run_id in dhandle.get_task_bold_run_ids(task)[subj]:
                        # run_id = 1
                        print "run{}".format(run_id)

                        # pad
                        if subj < 10:
                            subj_num = "sub00" + str(subj)
                        else:
                            subj_num = "sub0" + str(subj)

                        # determine mask
                        mask_fname = os.path.join(mkhd_data_path, subj_num,
                                                  "masks/orig/Group_Int_Mask.nii.gz")

                        # get model - this seems to be the problem - loads events from both tasks!
                        # run_events = dhandle.get_bold_run_model(
                        #     model, subj, run_id)

                        tmp_re = os.path.join(
                            mkhd_data_path, subj_num, "BOLD", "task00{}_run00{}".format(task, run_id), "attributes.txt")
                        run_events = SampleAttributes(tmp_re)
                        # print np.unique(run_events.targets)
                        # print len(run_events.targets)

                        # # get data
                        # run_data = dhandle.get_bold_run_dataset(subj,
                        #                                          task,
                        #                                          run_id,
                        #                                          chunks=run_id - 1,
                        #                                          mask=mask_fname)

                        # get behavior events
                        # run_data.sa['targets'] = events2sample_attr(run_events,
                        #                                              run_data.sa.time_coords,
                        #                                              noinfolabel='base')

                        tmp_rd = os.path.join(mkhd_data_path, subj_num, "BOLD",
                                              "task00{}_run00{}".format(task, run_id), "bold.nii.gz")
                        run_data = fmri_dataset(samples=tmp_rd, targets=run_events.targets,
                                                chunks=run_events.chunks, mask=mask_fname)
                        # print run_data.summary()

                        # clean, save
                        run_data = run_data[run_data.sa.targets != 'base']
                        run_data = run_data[run_data.sa.targets != 'num']
                        print run_data.summary()

                        subj_data.append(run_data)
                        del run_data

                    # # check that number of voxels (length) is equal across subjs
                    # #   this is a result of the mask
                    # print run_data.fa.values()

                    # append list, collapse across runs
                    group_data.append(vstack(subj_data, a=0))
                    del subj_data

                # save all data as one file, clean
                mvpa2.base.hdf5.h5save(out_file,
                                       group_data,
                                       mode='w',
                                       mkdir=True,
                                       compression="gzip")

                del group_data


# %%
# def func_train(hdf5_path, group_path):
"""
Step 2: Between, Within Classifiers
"""

# For testing
par_path = "/Users/nmuncy/Projects/learn_mvpa"
group_path = os.path.join(par_path, "grpAnalysis")
data_path = os.path.join(par_path, "mvpa")
hdf5_path = os.path.join(data_path, "hdf5")

# get data
fds_train = h5load(os.path.join(hdf5_path, "model1_data_Train.hdf5.gz"))

for i, sd in enumerate(fds_train):
    sd.sa['subject'] = np.repeat(i, len(sd))
nsubjs = len(fds_train)
ncats = len(fds_train[0].UT)
nruns = len(fds_train[0].UC)

# write out summary
h_out = """
    Number of Subjects: {}
    Number of Categories: {}
    Number of Runs: {}

    {}
""".format(nsubjs, ncats, nruns, fds_train[0].summary())
write_out = open(os.path.join(group_path, "mvpa_data_summary.txt"), "w")
write_out.write(h_out)
write_out.close()

# %%
# Set up feature selection
#   100 highest anova values
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
#   should errorfx be mean_match_accuracy? or lambda p, t: np.mean(p == t)

# cvte = CrossValidation(clf,
#                         NFoldPartitioner(),
#                         errorfx=lambda p, t: np.mean(p == t),
#                         enable_ca=['stats'])

cvws = CrossValidation(fsclf,
                       NFoldPartitioner(attr='chunks'),
                       errorfx=mean_match_accuracy)

cvbs = CrossValidation(fsclf,
                       NFoldPartitioner(attr='subject'),
                       errorfx=mean_match_accuracy)

# %%
# run classifier
#   TODO need to call tested samples the same

# get test data
fds_test = h5load(os.path.join(hdf5_path, "model1_data_Test.hdf5.gz"))

for i, sd in enumerate(fds_test):
    sd.sa['subject'] = np.repeat(i, len(sd))

# test
wsc_results = [cvws(sd) for sd in fds_test]
wsc_results = vstack(wsc_results)
fds_comb = vstack(fds_test)
bsc_results = cvbs(fds_comb)

# %%
"""
Step 3: Hyperalignment
"""
# classifier
cv = CrossValidation(clf,
                     NFoldPartitioner(attr='subject'),
                     errorfx=mean_match_accuracy)

bsc_hyper_results = []
for test_run in range(nruns):

    ds_train = [sd[sd.sa.chunks != test_run, :] for sd in fds_train]
    ds_test = [sd[sd.sa.chunks == test_run, :] for sd in fds_train]

    anova = OneWayAnova()
    fscores = [anova(sd) for sd in ds_train]
    featsels = [
        StaticFeatureSelection(fselector(fscore)) for fscore in fscores
    ]
    ds_train_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_train)]

    hyper = Hyperalignment()
    hypmaps = hyper(ds_train_fs)

    ds_test_fs = [fs.forward(sd) for fs, sd in zip(featsels, ds_test)]
    ds_hyper = [h.forward(sd) for h, sd in zip(hypmaps, ds_test_fs)]

    ds_hyper = vstack(ds_hyper)
    zscore(ds_hyper, chunks_attr='subject')
    res_cv = cv(ds_hyper)
    bsc_hyper_results.append(res_cv)

bsc_hyper_results = hstack(bsc_hyper_results)

# Compare Classifiers
h_out = """
    Classifier = LinearSVM
    Cross-validation = NFoldPartitioner (subject)

    Average Classification Accuracies
        Within-Subject: {} +/- {}
        Between-Subject: {} +/- {}
        Hyper Between-Subject: {} +/- {}
""".format(
    round(np.mean(wsc_results), 2),
    round(np.std(wsc_results) / np.sqrt(nsubjs - 1), 3),
    round(np.mean(bsc_results), 2),
    round(np.std(np.mean(bsc_results, axis=1)) / np.sqrt(nsubjs - 1), 3),
    round(np.mean(bsc_hyper_results), 2),
    round(
        np.std(np.mean(bsc_hyper_results, axis=1)) / np.sqrt(nsubjs - 1),
        3))

write_out = open(os.path.join(group_path, "mvpa_stats_classifiers.txt"),
                 "w")
write_out.write(h_out)
write_out.close()

# %%
# Similarity
anova = OneWayAnova()
fscores = [anova(sd) for sd in fds_train]
fscores = np.mean(np.asarray(vstack(fscores)), axis=0)

ds_fs = [sd[:, fselector(fscores)] for sd in fds_train]
hyper = Hyperalignment()
mappers = hyper(ds_fs)
ds_hyper = [m.forward(ds_) for m, ds_ in zip(mappers, ds_fs)]

sm_orig = [
    np.corrcoef(sd.get_mapped(mean_group_sample(['targets'])).samples)
    for sd in ds_fs
]
sm_orig_mean = np.mean(sm_orig, axis=0)

sm_hyper_mean = np.mean([
    np.corrcoef(sd.get_mapped(mean_group_sample(['targets'])).samples)
    for sd in ds_hyper
],
    axis=0)

ds_hyper = vstack(ds_hyper)
sm_hyper = np.corrcoef(ds_hyper.get_mapped(mean_group_sample(['targets'])))

ds_fs = vstack(ds_fs)
sm_anat = np.corrcoef(ds_fs.get_mapped(mean_group_sample(['targets'])))

intended_label_order = [0, 1]
labels = fds_train[0].UT
labels = labels[intended_label_order]

pl.figure(figsize=(12, 12))
# plot all three similarity structures
for i, sm_t in enumerate((
    (sm_orig_mean, "Average within-subject\nsimilarity"),
    (sm_anat, "Similarity of group average\ndata (anatomically aligned)"),
    (sm_hyper_mean,
        "Average within-subject\nsimilarity (hyperaligned data)"),
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
    pl.xticks(range(ncats),
              labels,
              size='small',
              stretch='ultra-condensed',
              rotation=45)
    pl.yticks(range(ncats),
              labels,
              size='small',
              stretch='ultra-condensed',
              rotation=45)
    pl.ylim(ylim)

pl.savefig(os.path.join(group_path, "mvpa_plot_classifiers.png"))


# %%
def main():

    # set up
    main_par_path = sys.argv[1]
    # main_par_path = "/Users/nmuncy/Projects/learn_mvpa"
    main_group_path = os.path.join(main_par_path, "grpAnalysis")
    main_data_path = os.path.join(main_par_path, "mvpa")
    main_hd5f_path = os.path.join(main_data_path, "hdf5")

    # get model dict
    with open(os.path.join(main_data_path, "model_dict.json")) as json_file:
        main_model_dict = json.load(json_file)

    # make hdf5 datasets
    if not os.path.exists(main_hd5f_path):
        os.makedirs(main_hd5f_path)
    func_make_hdf5(main_data_path, main_hd5f_path, main_model_dict)

    # run classifier
    # func_train(main_hd5f_path, main_group_path)


if __name__ == "__main__":
    main()

# %%
