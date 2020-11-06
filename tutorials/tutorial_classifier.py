# %%
from mvpa2.tutorial_suite import *

# %%
"""
Data In Shape
"""
# 2d set
ds = dataset_wizard(np.ones((5, 12)))
ds.shape
'mapper' in ds.a

# flat 3d
ds = dataset_wizard(np.ones((5, 4, 3)))
ds.shape
'mapper' in ds.a
print ds.a.mapper

# select features, mapper chain shows history
myfavs = [1, 2, 8, 10]
subds = ds[:, myfavs]
subds.shape
'mapper' in subds.a
print subds.a.mapper

# give data to existing mapper chain (subds.a.mapper)
fwdtest = np.arange(12).reshape(4, 3)
print fwdtest
fmapped = subds.a.mapper.forward1(fwdtest)
fmapped.shape
print fmapped


# %%
"""
Load real data (one subj/run)
    12 runs of viewing 8 image types + rest
    Target = Label = experiment condition / behavior
    Chunks = part of data
"""
data_path = "/home/nate/Projects/pymvpa/tutorial_data_4/data"

# attributes (behaviors) and chunks
#   this has 1 attribute per sample/volume
attr_fname = os.path.join(data_path, "sub001", "BOLD",
                          "task001_run001", "attributes.txt")
attr = SampleAttributes(attr_fname)
len(attr.targets)
print np.unique(attr.targets)

len(attr.chunks)
print np.unique(attr.chunks)

# mri data and attributes
bold_fname = os.path.join(data_path, "sub001", "BOLD",
                          "task001_run001", "bold.nii.gz")
mask_fname = os.path.join(data_path, "sub001", "masks", "orig", "vt.nii.gz")
fds = fmri_dataset(samples=bold_fname, targets=attr.targets,
                   chunks=attr.chunks, mask=mask_fname)
fds.shape
print fds.sa


# %%
"""
Load all runs of subj
    (should support multiple subj?)
    Supports attributes not aligned with volumes
    Note dir organization
"""
# load data
dhandle = OpenFMRIDataset(data_path)
dhandle.get_subj_ids()
dhandle.get_task_descriptions()

# get timing info from subj/model files (timing files)
model = 1
subj = 1
run = 1
events = dhandle.get_bold_run_model(model, subj, run)
for ev in events[:2]:
    print ev

# integrate timing (events * sample acquisition)
targets = events2sample_attr(
    events, fds.sa.time_coords, noinfolabel='rest', onset_shift=0.0)
print np.unique([attr.targets[i] == t for i, t in enumerate(targets)])
print np.unique(attr.targets)
print len(fds), len(targets)

# access fmri data
task = 1
fds = dhandle.get_bold_run_dataset(subj, task, run, mask=mask_fname)
print fds


# %%
"""
Multi-session Data
"""
# combine all runs into run_datasets
task = 1
model = 1
subj = 1
run_datasets = []

for run_id in dhandle.get_task_bold_run_ids(task)[subj]:

    # design info
    run_events = dhandle.get_bold_run_model(model, subj, run_id)

    # fmri data, mask
    run_ds = dhandle.get_bold_run_dataset(
        subj, task, run_id, chunks=run_id - 1, mask=mask_fname)

    # convert event to sample attribute, assign as target
    run_ds.sa['targets'] = events2sample_attr(
        run_events, run_ds.sa.time_coords, noinfolabel='rest')

    # write list
    run_datasets.append(run_ds)

# merge datasets, a=0 means attributes for first
#   set should be used for all data
fds = vstack(run_datasets, a=0)
print fds.summary()


# %%
"""
Basic Preprocessing
"""
detrender = PolyDetrendMapper(polyord=1, chunks_attr='chunks')
detrended_fds = fds.get_mapped(detrender)
print detrended_fds.a.mapper

zscorer = ZScoreMapper(param_est=('targets', ['rest']))

zscore(detrended_fds, param_est=('targets', ['rest']))
fds = detrended_fds
print fds.a.mapper

fds = fds[fds.sa.targets != 'rest']
print fds.shape


# %%
"""
Computing Patterns of Activation
    This does mean activity per target
"""
rnames = {0: 'even', 1: 'odd'}
fds.sa['runtype'] = [rnames[c % 2] for c in fds.sa.chunks]

averager = mean_group_sample(['targets', 'runtype'])
type(averager)
fds = fds.get_mapped(averager)
fds.shape
print fds.sa.targets
print fds.sa.chunks


# %%
"""
There and back again
    Data from first section
    Move in both directions of mapper
"""
print ds
print ds.a.mapper

# move backwards
orig_data = ds.a.mapper.reverse(ds.samples)
orig_data.shape

# o/syntax
orig_data = ds.O
orig_data.shape

# reversing along chain
print subds
print subds.a.mapper
subds.nfeatures
revtest = np.arange(subds.nfeatures) + 10
print revtest
rmapped = subds.a.mapper.reverse1(revtest)
rmapped.shape
print rmapped

# reverse fMRI
print fds.a.mapper
fds.nfeatures
revtest = np.arange(100, 100 + fds.nfeatures)
rmapped = fds.a.mapper.reverse1(revtest)
rmapped.shape

rmapped_partial = fds.a.mapper[:2].reverse1(revtest)
(rmapped == rmapped_partial).all()

# %%
# back to nifti
'imghdr' in fds.a
nimg = map2nifti(fds, revtest)
# nimg.to_filename('mytest.nii.gz')


# %%
"""
Classifiers
    kNN = nearest-neighbor classifier with correlation as
        a distance measure
        i.e. similarity of sample to training
"""
ds = fds

# set classifier (clf) as kNN, train
clf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
clf.train(ds)

# compare predictions (to training data)
predictions = clf.predict(ds.samples)
np.mean(predictions == ds.sa.targets)

# split data to training, sample sets
print ds.sa.runtype
ds_split1 = ds[ds.sa.runtype == 'odd']
len(ds_split1)
ds_split2 = ds[ds.sa.runtype == 'even']
len(ds_split2)

# rather than repeat above, call classifier
#   compute mismatch bx prediction and stored target values
#   train/test on both halves
clf.set_postproc(BinaryFxNode(mean_mismatch_error, 'targets'))

clf.train(ds_split2)
err = clf(ds_split1)
print np.asscalar(err)

clf.train(ds_split1)
err = clf(ds_split2)
print np.asscalar(err)

# %%
"""
Cross-validation
    Half split is one method
        HalfPartitioner can split data for us
    Leave-one-out is another method

"""
clf.set_postproc(None)
hpart = HalfPartitioner(attr='runtype')
cv = CrossValidation(clf, hpart)
cv_results = cv(ds)

np.mean(cv_results)
len(cv_results)
cv_results.samples

# change classifier type
clf = LinearCSVMC()
cvte = CrossValidation(clf, HalfPartitioner(attr='runtype'))
cv_results = cvte(ds)
np.mean(cv_results)

# get accuracy (rather than error)
cvte = CrossValidation(clf, HalfPartitioner(
    attr='runtype'), errorfx=lambda p, t: np.mean(p == t))
cv_results = cvte(ds)
np.mean(cv_results)

# %%
# --- Leave one run out
#   Get mean sample per run (not per target)

# get, preproc new data
ds = vstack(run_datasets, a=0)
poly_detrend(ds, polyord=1, chunks_attr='chunks')
zscore(ds, param_est=('targets', ['rest']))
ds = ds[ds.sa.targets != 'rest']

run_averager = mean_group_sample(['targets', 'chunks'])
ds = ds.get_mapped(run_averager)
ds.shape

# different cross-validation approach, get accuracy
cvte = CrossValidation(clf, NFoldPartitioner(),
                       errorfx=lambda p, t: np.mean(p == t))
cv_results = cvte(ds)
np.mean(cv_results)

type(cv_results)
print cv_results.samples


# %%
"""
Confusion matrices
    Get accuracy for e/target, not average across
        types in each run
    confusion matrix stored on conditional attributes (ca)
        collection
"""
cvte = CrossValidation(clf, NFoldPartitioner(
), errorfx=lambda p, t: np.mean(p == t), enable_ca=['stats'])
cv_results = cvte(ds)
print cvte.ca.stats.as_string(description=True)
print cvte.ca.stats.matrix
