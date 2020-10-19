import shutil
import tempfile
import os
from glob import glob
import numpy as np
from mvpa2.tutorial_suite import *

"""
Datasets
"""
# first axis = row = sample = fMRI volumes
# second axis = column = feature = voxels
# (sample, feature)
data = [[1, 1, -1], [2, 0, 0], [3, 1, 1], [4, 0, -1]]
ds = Dataset(data)
ds.shape
len(ds)
ds.nfeatures
ds.samples

one_d = [0, 1, 2, 3]
one_ds = Dataset(one_d)
one_ds.shape

# for multi-dimensional data, only 2nd value is feature
#   (x, feature, y, z)
#   (3, 4, 2, 3) = three samples of four features,
#       e/sample is 2x3 matrix
m_ds = Dataset(np.random.random((3, 4, 2, 3)))
m_ds.shape
m_ds.nfeatures


"""
Attributes
"""
# --- Samples
#   accessible via .sa
#   similar to dict
#   one value per sample
ds.sa['some_attr'] = [0., 1, 1, 3]
ds.sa.keys()

# samples stored as "collectable", not text
#   allows for quick access of values
type(ds.sa['some_attr'])
ds.sa['some_attr'].value
ds.sa['some_attr'].unique
ds.sa.some_attr

# forces format
ds.sa['invalid'] = 4
ds.sa['invalid'] = [1, 2, 3, 4, 5, 6]

# can store different data types in same dataframe
ds.sa['literal'] = ['one', 'two', 'three', 'four']
sorted(ds.sa.keys())
for attr in ds.sa:
    print "%s: %s" % (attr, ds.sa[attr].value.dtype.name)


# --- Features
#   accessible via .fa
#   one value per feature
ds.nfeatures
ds.fa['my_fav'] = [0, 1, 0]
ds.fa['responsible'] = ['me', 'you', 'nobody']
sorted(ds.fa.keys())


# --- Attributes
#   about whole dataset
#   accessible via .a
#   no constrains - anything can be stored
ds.a['pointless'] = glob("*/*")
'setup.py' in ds.a.pointless


"""
Slicing
"""
ds.samples

# multiple ways of slicing (every other)
ds[::2].samples

mask = np.array([True, False, True, False])
ds[mask].samples

ds[[0, 2]].samples

# features on o/axis
#   get only 2 features of all samples
ds[:, [1, 2]].samples

# slice for sample & feature
subds = ds[[0, 1], [0, 2]]
subds.samples

# what about attributes after slicing for
#   samples/features?
print ds.sa.some_attr
print ds.fa.responsible
print subds.sa.some_attr
print subds.fa.responsible

# selecting off of .sa or .fa
subds = ds[ds.sa.some_attr == 1., ds.fa.responsible == 'me']
print subds.shape

subds = ds[{'some_attr': [1., 0.], 'literal': ['two']},
           {'responsible': ['me', 'you']}]
print subds.sa.some_attr, subds.sa.literal, subds.fa.responsible


"""
Load fMRI data
"""
tutorial_data_path = "/home/nate/Projects/pymvpa/tutorial_data_4"
bold_fname = os.path.join(tutorial_data_path, "data",
                          "sub001", "BOLD", "task001_run001", "bold.nii.gz")

# samples = volumes, features = voxels (collapsed to 1D)
ds = fmri_dataset(bold_fname)
len(ds)
ds.nfeatures
ds.shape

# masking
mask_fname = os.path.join(tutorial_data_path, "data",
                          "sub001", "masks", "orig", "vt.nii.gz")
ds = fmri_dataset(bold_fname, mask=mask_fname)
len(ds)
ds.nfeatures

# sample/volume info
ds.sa.time_indices[:5]
ds.sa.time_coords[:5]

# voxel info
ds.fa.voxel_indices[:5]
ds.a.voxel_eldim
ds.a.voxel_dim
'imghdr' in ds.a

# mapper
print ds.a.mapper

# get some info only
stripped = ds.copy(deep=False, sa=['time_coords'], fa=[], a=[])
print stripped


"""
Intermediate storage
"""
tempdir = tempfile.mkdtemp()
ds.save(os.path.join(tempdir, 'mydataset.hdf5'))

ds.save(os.path.join(tempdir, 'mydataset.gzipped.hdf5'), compression=9)
h5save(os.path.join(tempdir, 'mydataset.gzipped.hdf5'), ds, compression=9)

loaded = h5load(os.path.join(tempdir, 'mydataset.hdf5'))
np.all(ds.samples == loaded.samples)
shutil.rmtree(tempdir, ignore_errors=True)
