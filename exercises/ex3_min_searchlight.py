import numpy as np

from mvpa2.generators.partition import OddEvenPartitioner
from mvpa2.clfs.svm import LinearCSVMC
from mvpa2.measures.base import CrossValidation
from mvpa2.measures.searchlight import sphere_searchlight
from mvpa2.testing.datasets import datasets
from mvpa2.mappers.fx import mean_sample

dataset = datasets['3dlarge']

# setup measure to be computed in each sphere (cross-validated
# generalization error on odd/even splits)
cv = CrossValidation(LinearCSVMC(), OddEvenPartitioner())

# setup searchlight with 2 voxels radius and measure configured above
sl = sphere_searchlight(cv, radius=2, space='myspace',
                        postproc=mean_sample())

# run searchlight on dataset
sl_map = sl(dataset)

print 'Best performing sphere error:', np.min(sl_map.samples)
