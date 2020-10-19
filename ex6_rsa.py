# %%
from nipy.labs.viz_tools.activation_maps import plot_map
import nibabel as nb
from mvpa2.mappers.shape import TransposeMapper
from mvpa2.base.learner import ChainLearner
from scipy.spatial.distance import squareform
from mvpa2.measures.searchlight import sphere_searchlight
from mvpa2.measures import rsa
from mvpa2.mappers.fx import mean_group_sample
from mvpa2.mappers.zscore import zscore
from mvpa2.mappers.detrend import poly_detrend
import numpy as np
import pylab as pl
from os.path import join as pjoin
from mvpa2.suite import *
from mvpa2 import cfg
from mvpa2.datasets.sources.native import load_tutorial_data

datapath = "/home/nate/Projects/pymvpa/tutorial_data_4/data"
mask_fname = os.path.join(datapath, "sub001", "masks", "orig", "hoc.nii.gz")

# ds = load_tutorial_data(roi=(15, 16, 23, 24, 36, 38, 39, 40, 48))
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

ds = vstack(run_datasets, a=0)


poly_detrend(ds, polyord=1, chunks_attr='chunks')
# z-scoring with respect to the 'rest' condition
zscore(ds, chunks_attr='chunks', param_est=('targets', 'rest'))
# now remove 'rest' samples
ds = ds[ds.sa.targets != 'rest']

# %%


def plot_mtx(mtx, labels, title):
    pl.figure()
    pl.imshow(mtx, interpolation='nearest')
    pl.xticks(range(len(mtx)), labels, rotation=-45)
    pl.yticks(range(len(mtx)), labels)
    pl.title(title)
    pl.clim((0, 2))
    pl.colorbar()


# compute a dataset with the mean samples for all conditions
mtgs = mean_group_sample(['targets'])
mtds = mtgs(ds)

# basic ROI RSA -- dissimilarity matrix for the entire ROI
dsm = rsa.PDist(square=True)
res = dsm(mtds)
plot_mtx(res, mtds.sa.targets, 'ROI pattern correlation distances')

# same as above, but done in a searchlight fashion
dsm = rsa.PDist(square=False)
sl = sphere_searchlight(dsm, 2)
slres = sl(mtds)

# score each searchlight sphere result wrt global pattern dissimilarity
distinctiveness = np.sum(np.abs(slres), axis=0)
print 'Most dissimilar patterns around', \
    mtds.fa.voxel_indices[distinctiveness.argmax()]
# take a look at the this dissimilarity structure
plot_mtx(squareform(slres.samples[:, distinctiveness.argmax()]),
         mtds.sa.targets,
         'Maximum distinctive searchlight pattern correlation distances')

# %%
# more interesting: let's look at the stability of similarity structures
# across experiment runs
# mean condition samples, as before, but now individually for each run
mtcgs = mean_group_sample(['targets', 'chunks'])
mtcds = mtcgs(ds)

# searchlight consistency measure -- how correlated are the structures
# across runs
dscm = rsa.PDistConsistency()
sl_cons = sphere_searchlight(dscm, 2)
slres_cons = sl_cons(mtcds)

# mean correlation
mean_consistency = np.mean(slres_cons, axis=0)
print 'Most stable dissimilarity patterns around', \
    mtds.fa.voxel_indices[mean_consistency.argmax()]
# Look at this pattern
plot_mtx(squareform(slres.samples[:, mean_consistency.argmax()]),
         mtds.sa.targets,
         'Most consistent searchlight pattern correlation distances')

# %%
# let's see where in the brain we find dissimilarity structures that are
# similar to our most stable one
tdsm = rsa.PDistTargetSimilarity(
    slres.samples[:, mean_consistency.argmax()])
# using a searchlight
sl_tdsm = sphere_searchlight(ChainLearner([tdsm, TransposeMapper()]), 2)
slres_tdsm = sl_tdsm(mtds)

# %%
# plot the spatial distribution using NiPy
vol = ds.a.mapper.reverse1(slres_tdsm.samples[0])
anat = nb.load(pjoin(datapath, 'sub001', 'anatomy', 'highres001.nii.gz'))

pl.figure(figsize=(15, 4))
sp = pl.subplot(121)
pl.title('Distribution of target similarity structure correlation')
slices = plot_map(
    vol,
    ds.a.imgaffine,
    cut_coords=np.array((12, -42, -20)),
    threshold=.5,
    cmap="bwr",
    vmin=0,
    vmax=1.,
    axes=sp,
    anat=anat.get_data(),
    anat_affine=anat.affine,
)
img = pl.gca().get_images()[1]
cax = pl.axes([.05, .05, .05, .9])
pl.colorbar(img, cax=cax)

sp = pl.subplot(122)
pl.hist(slres_tdsm.samples[0],
        normed=False,
        bins=30,
        color='0.6')
pl.ylabel("Number of voxels")
pl.xlabel("Target similarity structure correlation")

# %%
