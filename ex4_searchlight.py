# %%

from mvpa2.suite import *

# enable debug output for searchlight call
if __debug__:
    debug.active += ["SLC"]

# get data

# dataset = load_tutorial_data(
#     roi='brain',
#     add_fa={'vt_thr_glm': os.path.join(datapath, 'sub001', 'masks',
#                                        'orig', 'vt.nii.gz')})

datapath = "/home/nate/Projects/pymvpa/tutorial_data_4/data"
mask_fname = os.path.join(datapath, "sub001", "masks", "orig", "vt.nii.gz")
dhandle = OpenFMRIDataset(datapath)

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

dataset = vstack(run_datasets, a=0)
# print dataset.summary()
poly_detrend(dataset, polyord=1, chunks_attr='chunks')
dataset = dataset[np.array([l in ['rest', 'house', 'scrambledpix']
                            for l in dataset.targets], dtype='bool')]
zscore(dataset, chunks_attr='chunks', param_est=(
    'targets', ['rest']), dtype='float32')

dataset = dataset[dataset.sa.targets != 'rest']

# %%
# choose classifier
clf = LinearNuSVMC()

# setup measure to be computed by Searchlight
# cross-validated mean transfer using an N-fold dataset splitter
cv = CrossValidation(clf, NFoldPartitioner())

# center_ids = dataset.fa.vt_thr_glm.nonzero()[0]

# for radius in [0, 1, 3]:
#     # tell which one we are doing
#     print "Running searchlight with radius: %i ..." % (radius)

#     sl = sphere_searchlight(cv, radius=radius, space='voxel_indices',
#                             postproc=mean_sample())

#     ds = dataset.copy(deep=False,
#                         sa=['targets', 'chunks'],
#                         fa=['voxel_indices'],
#                         a=['mapper'])

#     sl_map = sl(ds)
#     sl_map.samples *= -1
#     sl_map.samples += 1

radius = 3
sl = sphere_searchlight(cv, radius=radius, space='voxel_indices',
                        postproc=mean_sample())

ds = dataset.copy(deep=False,
                  sa=['targets', 'chunks'],
                  fa=['voxel_indices'],
                  a=['mapper'])

sl_map = sl(ds)
sl_map.samples *= -1
sl_map.samples += 1

# %%
niftiresults = map2nifti(sl_map, imghdr=dataset.a.imghdr)

# setup plotting parameters (not essential for the analysis itself)
plot_args = {
    'background': os.path.join(datapath, 'sub001', 'anatomy', 'highres001.nii.gz'),
    'background_mask': os.path.join(datapath, 'sub001', 'masks', 'orig', 'brain.nii.gz'),
    'overlay_mask': os.path.join(datapath, 'sub001', 'masks', 'orig', 'vt.nii.gz'),
    'do_stretch_colors': False,
    'cmap_bg': 'gray',
    'cmap_overlay': 'autumn',  # YlOrRd_r # pl.cm.autumn
    'interactive': cfg.getboolean('examples', 'interactive', True),
}

fig = pl.figure(figsize=(12, 12), facecolor='white')
subfig = plot_lightbox(overlay=niftiresults,
                       vlim=(0.5, None), slices=range(23, 31),
                       fig=fig, **plot_args)
pl.title('Accuracy distribution for radius %i' % radius)
