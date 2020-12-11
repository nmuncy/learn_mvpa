"""
Notes

Written for Python 2.7 <-------------!!!

Also written for local environment
"""


# %%
from mvpa2.suite import *
import os
import gc
import sys
import json


"""
Step 0: Set up
"""
par_path = "/Users/nmuncy/Projects/learn_mvpa"
group_path = os.path.join(par_path, "grpAnalysis")
hdf5_path = os.path.join(par_path, "mvpa/hdf5")

if not os.path.exists(group_path):
    os.makedirs(group_path)


# %%
"""
Step 1: Train Between, Within Classifiers
"""
# get data
fds_train = h5load(os.path.join(hdf5_path, "model1_data_Train.hdf5.gz"))

for i, sd in enumerate(fds_train):
    sd.sa["subject"] = np.repeat(i, len(sd))
train_nsubj = len(fds_train)
train_ncats = len(fds_train[0].UT)
train_nruns = len(fds_train[0].UC)

# write out summary
h_out = """
    Number of Subjects: {}
    Number of Categories: {}
    Number of Runs: {}

    {}
""".format(
    train_nsubj, train_ncats, train_nruns, fds_train[0].summary()
)
write_out = open(os.path.join(group_path, "mvpa_train_summary.txt"), "w")
write_out.write(h_out)
write_out.close()

# %%
# Set up feature selection
#   100 highest anova values
clf = LinearCSVMC()
nf = 100
fselector = FixedNElementTailSelector(
    nf, tail="upper", mode="select", sort=False)

sbfs = SensitivityBasedFeatureSelection(
    OneWayAnova(), fselector, enable_ca=["sensitivities"]
)

fsclf = FeatureSelectionClassifier(clf, sbfs)

# Classify within, between
#   should errorfx be mean_match_accuracy? or lambda p, t: np.mean(p == t)
# cvte = CrossValidation(clf,
#                         NFoldPartitioner(),
#                         errorfx=lambda p, t: np.mean(p == t),
#                         enable_ca=['stats'])


# Don't cross validate
# # within
# cvws = CrossValidation(
#     fsclf, NFoldPartitioner(attr="chunks"), errorfx=mean_match_accuracy
# )

# # between
# cvbs = CrossValidation(
#     fsclf, NFoldPartitioner(attr="subject"), errorfx=mean_match_accuracy
# )

fsclf.train(fds_train[0])


# %%
"""
Step 2: Test classifiers
"""
# get test data
fds_test = h5load(os.path.join(hdf5_path, "model1_data_Test.hdf5.gz"))

for i, sd in enumerate(fds_test):
    sd.sa["subject"] = np.repeat(i, len(sd))
test_nsubj = len(fds_test)
test_ncats = len(fds_test[0].UT)
test_nruns = len(fds_test[0].UC)

# write out summary
h_out = """
    Number of Subjects: {}
    Number of Categories: {}
    Number of Runs: {}

    {}
""".format(
    test_nsubj, test_ncats, test_nruns, fds_test[0].summary()
)
write_out = open(os.path.join(group_path, "mvpa_test_summary.txt"), "w")
write_out.write(h_out)
write_out.close()


# test - does this work?
# wsc_results = [cvws(sd) for sd in fds_test]
wsc_results = [fsclf(x) for x in fds_test]
wsc_results = vstack(wsc_results)
fds_comb = vstack(fds_test)

# bsc_results = cvbs(fds_comb)
bsc_results = fsclf(fds_comb)


# %%
"""
Step 3: Train, Test Hyperalignment
"""
# classifier
cv = CrossValidation(clf, NFoldPartitioner(
    attr="subject"), errorfx=mean_match_accuracy)

bsc_hyper_results = []
for test_run in range(nruns):

    ds_train = [sd[sd.sa.chunks != test_run, :] for sd in fds_train]
    ds_test = [sd[sd.sa.chunks == test_run, :] for sd in fds_train]

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
    zscore(ds_hyper, chunks_attr="subject")
    res_cv = cv(ds_hyper)
    bsc_hyper_results.append(res_cv)

bsc_hyper_results = hstack(bsc_hyper_results)


# %%
"""
Step 4: Compare Classifiers
"""
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
    round(np.std(np.mean(bsc_hyper_results, axis=1)) / np.sqrt(nsubjs - 1), 3),
)

write_out = open(os.path.join(group_path, "mvpa_stats_classifiers.txt"), "w")
write_out.write(h_out)
write_out.close()


# %%
"""
Step 5: Show Similarity
"""
anova = OneWayAnova()
fscores = [anova(sd) for sd in fds_train]
fscores = np.mean(np.asarray(vstack(fscores)), axis=0)

ds_fs = [sd[:, fselector(fscores)] for sd in fds_train]
hyper = Hyperalignment()
mappers = hyper(ds_fs)
ds_hyper = [m.forward(ds_) for m, ds_ in zip(mappers, ds_fs)]

sm_orig = [
    np.corrcoef(sd.get_mapped(mean_group_sample(["targets"])).samples) for sd in ds_fs
]
sm_orig_mean = np.mean(sm_orig, axis=0)

sm_hyper_mean = np.mean(
    [
        np.corrcoef(sd.get_mapped(mean_group_sample(["targets"])).samples)
        for sd in ds_hyper
    ],
    axis=0,
)

ds_hyper = vstack(ds_hyper)
sm_hyper = np.corrcoef(ds_hyper.get_mapped(mean_group_sample(["targets"])))

ds_fs = vstack(ds_fs)
sm_anat = np.corrcoef(ds_fs.get_mapped(mean_group_sample(["targets"])))

intended_label_order = [0, 1]
labels = fds_train[0].UT
labels = labels[intended_label_order]

pl.figure(figsize=(12, 12))
# plot all three similarity structures
for i, sm_t in enumerate(
    (
        (sm_orig_mean, "Average within-subject\nsimilarity"),
        (sm_anat, "Similarity of group average\ndata (anatomically aligned)"),
        (sm_hyper_mean, "Average within-subject\nsimilarity (hyperaligned data)"),
        (sm_hyper, "Similarity of group average\ndata (hyperaligned)"),
    )
):
    sm, title = sm_t
    # reorder matrix columns to match label order
    sm = sm[intended_label_order][:, intended_label_order]
    pl.subplot(2, 2, i + 1)
    pl.imshow(sm, vmin=-1.0, vmax=1.0, interpolation="nearest")
    pl.colorbar(shrink=0.4, ticks=[-1, 0, 1])
    pl.title(title, size=12)
    ylim = pl.ylim()
    pl.xticks(
        range(ncats), labels, size="small", stretch="ultra-condensed", rotation=45
    )
    pl.yticks(
        range(ncats), labels, size="small", stretch="ultra-condensed", rotation=45
    )
    pl.ylim(ylim)

pl.savefig(os.path.join(group_path, "mvpa_plot_classifiers.png"))
