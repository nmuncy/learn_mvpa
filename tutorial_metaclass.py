# %%
from mvpa2.tutorial_suite import *


"""
Meta-Classifier
    Use different classifiers on same data
"""

# get data
data_path = "/home/nate/Projects/pymvpa/tutorial_data_4/data"
mask_fname = os.path.join(data_path, "sub001", "masks", "orig", "vt.nii.gz")

dhandle = OpenFMRIDataset(data_path)
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
ds = vstack(run_datasets, a=0)

poly_detrend(ds, polyord=1, chunks_attr='chunks')
zscore(ds, param_est=('targets', ['rest']))
ds = ds[ds.sa.targets != 'rest']

run_averager = mean_group_sample(['targets', 'chunks'])
ds = ds.get_mapped(run_averager)
ds.shape

# %%
# classify not voxel intensity, but samples
#   in space spanned by single vectors
baseclf = LinearCSVMC()
metaclf = MappedClassifier(baseclf, SVDMapper())
cvte = CrossValidation(metaclf, NFoldPartitioner())
cv_results = cvte(ds)
print np.mean(cv_results)

# specify feature (see how errors increase)
mapper = ChainMapper([SVDMapper(), StaticFeatureSelection(slice(None, 2))])
metaclf = MappedClassifier(baseclf, mapper)
cvte = CrossValidation(metaclf, NFoldPartitioner())
cv_results = cvte(ds)
svm_err = np.mean(cv_results)
print round(svm_err, 2)

# try different classifier (doesn't help)
baseclf = kNN(k=1, dfx=one_minus_correlation, voting='majority')
mapper = ChainMapper([SVDMapper(), StaticFeatureSelection(slice(None, 2))])
metaclf = MappedClassifier(baseclf, mapper)
cvte = CrossValidation(metaclf, NFoldPartitioner())
cv_results = cvte(ds)
np.mean(cv_results) < svm_err

# add confusion matrix (of kNN)
cvte = CrossValidation(metaclf, NFoldPartitioner(), enable_ca=['stats'])
cv_results = cvte(ds)
print cvte.ca.stats.as_string(description=True)
print cvte.ca.stats.matrix

cv_acc = cv_results.samples
hist(cv_acc, bins=np.linspace(0, 1, 18))
