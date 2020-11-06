from mvpa2.suite import *

# load PyMVPA example dataset with literal labels
dataset = load_example_fmri_dataset(literal=True)

poly_detrend(dataset, polyord=1, chunks_attr='chunks')
dataset = dataset[np.array([l in ['face', 'house']
                            for l in dataset.sa.targets], dtype='bool')]
cv = CrossValidation(SMLR(), OddEvenPartitioner(), errorfx=mean_mismatch_error)
error = cv(dataset)
print "Error for %i-fold cross-validation on %i-class problem: %f" % (
    len(dataset.UC), len(dataset.UT), np.mean(error))
