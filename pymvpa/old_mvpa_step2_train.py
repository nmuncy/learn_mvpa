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
    run_datasets = []

    for run_id in dhandle.get_task_bold_run_ids(task)[subj]:

        # design info
        run_events = dhandle.get_bold_run_model(model, subj, run_id)

        # fmri data, mask
        run_ds = dhandle.get_bold_run_dataset(subj,
                                              task,
                                              run_id,
                                              chunks=run_id - 1,
                                              mask=mask_subj)

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
def func_make_group(data_path, hdf5_path):
    """
    MVPA with all data

    Step 1: Load data
    """
    # # For testing
    # par_path = "/Users/nmuncy/Projects/learn_mvpa"
    # data_path = os.path.join(h_par_path, "mvpa")
    # hdf5_path = os.path.join(data_path, "hdf5")

    # load dataset
    dhandle = mvpa2.datasets.sources.OpenFMRIDataset(data_path)
    dhandle.get_subj_ids()
    dhandle.get_task_descriptions()

    task = 1
    model = 1
    data_list = []

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

            # clean, save
            run_ds = run_ds[run_ds.sa.targets != 'base']
            run_ds = run_ds[run_ds.sa.targets != 'num']
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


def func_train(hdf5_path, group_path):
    # %%
    """
    Step 2: Between, Within Classifiers
    """

    # For testing
    # par_path = "/Users/nmuncy/Projects/learn_mvpa"
    # group_path = os.path.join(par_path, "grpAnalysis")
    # data_path = os.path.join(par_path, "mvpa")
    # hdf5_path = os.path.join(data_path, "hdf5")

    # get data
    fds_all = h5load(os.path.join(hdf5_path, "data_all.hdf5.gz"))

    for i, sd in enumerate(fds_all):
        sd.sa['subject'] = np.repeat(i, len(sd))
    nsubjs = len(fds_all)
    ncats = len(fds_all[0].UT)
    nruns = len(fds_all[0].UC)

    # write out summary
    h_out = """
        Number of Subjects: {}
        Number of Categories: {}
        Number of Runs: {}

        {}
    """.format(nsubjs, ncats, nruns, fds_all[0].summary())
    write_out = open(os.path.join(group_path, "mvpa_data_summary.txt"), "w")
    write_out.write(h_out)
    write_out.close()

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

    # run classifier
    wsc_results = [cvws(sd) for sd in fds_all]
    wsc_results = vstack(wsc_results)
    fds_comb = vstack(fds_all)
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

        ds_train = [sd[sd.sa.chunks != test_run, :] for sd in fds_all]
        ds_test = [sd[sd.sa.chunks == test_run, :] for sd in fds_all]

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
    fscores = [anova(sd) for sd in fds_all]
    fscores = np.mean(np.asarray(vstack(fscores)), axis=0)

    ds_fs = [sd[:, fselector(fscores)] for sd in fds_all]
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
    labels = fds_all[0].UT
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
    h_par_path = sys.argv[1]
    # h_par_path = "/Users/nmuncy/Projects/learn_mvpa"
    h_group_path = os.path.join(h_par_path, "grpAnalysis")
    h_data_path = os.path.join(h_par_path, "mvpa")
    h_hdf5_path = os.path.join(h_data_path, "hdf5")

    if not os.path.exists(h_hdf5_path):
        os.makedirs(h_hdf5_path)

    # load all subject data
    if not os.path.exists(os.path.join(h_hdf5_path, "data_all.hdf5.gz")):
        func_make_group(h_data_path, h_hdf5_path)

    # run classifier
    func_train(h_hdf5_path, h_group_path)


if __name__ == "__main__":
    main()

# %%
