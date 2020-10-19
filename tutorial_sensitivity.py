from mvpa2.tutorial_suite import *

"""
Sensitivity
    Train on only relevant voxels
    So, identify voxels with anova,
        then train on those.
    Train/Test must be independent
"""
# Get data
data_path = "/home/nate/Projects/pymvpa/tutorial_data_4/data"
mask_fname = os.path.join(data_path, "sub001", "masks", "orig", "brain.nii.gz")

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
print ds.shape

poly_detrend(ds, polyord=1, chunks_attr='chunks')
zscore(ds, param_est=('targets', ['rest']))
ds = ds[ds.sa.targets != 'rest']
run_averager = mean_group_sample(['targets', 'chunks'])
ds = ds.get_mapped(run_averager)
print ds.shape

# train classifier
clf = LinearCSVMC()
cvte = CrossValidation(clf, NFoldPartitioner(), enable_ca=['stats'])
results = cvte(ds)

# check classifier (wasn't great)
print np.round(cvte.ca.stats.stats['ACC%'], 1)
print cvte.ca.stats.matrix

# focus in on a few features only
#   i.e. find voxels w/diff signal for A vs B
#   first, get f-stats, find 500 largest f-values
fsel = SensitivityBasedFeatureSelection(
    OneWayAnova(), FixedNElementTailSelector(500, mode='select', tail='upper'))

fsel.train(ds)
ds_p = fsel(ds)
print ds_p.shape

# now run classifier on reduced dataset
#   stronger diagonal results
#   from double dipping
results = cvte(ds_p)
print np.round(cvte.ca.stats.stats['ACC%'], 1)
print cvte.ca.stats.matrix

# subset features
bin_demo = ds[np.array([i in ['bottle', 'shoe'] for i in ds.sa.targets])]
results = cvte(bin_demo)
print np.round(cvte.ca.stats.stats['ACC%'], 1)

fsel.train(bin_demo)
bin_demo_p = fsel(bin_demo)
results = cvte(bin_demo_p)
print cvte.ca.stats.stats["ACC%"]

# %%
# Do it correctly
#   i.e. no double-dipping
#   Use function to identify features
#       clf = classification function
#       fsel = f-stat selector
fclf = FeatureSelectionClassifier(clf, fsel)
cvte = CrossValidation(fclf, NFoldPartitioner(), enable_ca=['stats'])

# test on subset data
results = cvte(bin_demo)
print np.round(cvte.ca.stats.stats['ACC%'], 1)

# run on full data
results = cvte(ds)
print np.round(cvte.ca.stats.stats['ACC%'], 1)
print cvte.ca.stats.matrix

# %%
# Classify more features - top 5%
fsel = SensitivityBasedFeatureSelection(
    OneWayAnova(), FractionTailSelector(0.05, mode='select', tail='upper'))
fclf = FeatureSelectionClassifier(clf, fsel)
cvte = CrossValidation(fclf, NFoldPartitioner(), enable_ca=['stats'])
results = cvte(ds)
print np.round(cvte.ca.stats.stats['ACC%'], 1)

# get the weights
#   gives 28 sensitivty maps (all binary problems of 8 features)
sensana = fclf.get_sensitivity_analyzer()
type(sensana)
sens = sensana(ds)
type(sens)
print sens.shape

# combine
sens_comb = sens.get_mapped(maxofabs_sample())

# %%
# repeat, with cross-validation
sensana = fclf.get_sensitivity_analyzer(postproc=maxofabs_sample())

cv_sensana = RepeatedMeasure(sensana, ChainNode(
    (NFoldPartitioner(), Splitter('partitions', attr_values=(1,)))))
sens = cv_sensana(ds)
print sens.shape

ov = MapOverlap()
overlap_fraction = ov(sens.samples > 0)

# get matrix
sclf = SplitClassifier(fclf, enable_ca=['stats'])
cv_sensana = sclf.get_sensitivity_analyzer()
sens = cv_sensana(ds)
print sens.shape
print cv_sensana.clf.ca.stats.matrix
