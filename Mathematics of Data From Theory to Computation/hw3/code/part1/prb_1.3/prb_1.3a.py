import scipy.io
import numpy as np
from scipy.sparse import csr_matrix
from proj_L1 import proj_L1
from time import time


def proj_nuc(Z, kappa):
    #proj_nuc: This function implements the projection onto nuclear norm ball.

    # Implement projection operator here!

    assert kappa > 0, f"The radius should be positive, but the current radius is {kappa}"

    U, sigma, VT = np.linalg.svd(Z, full_matrices=False)  # SVD
    nuclear_norm = np.sum(sigma)  # nuclear norm
    # check if Z itself is a solution
    if nuclear_norm <= kappa:
        return Z
    return U @ np.diag(proj_L1(sigma)) @ VT  # return projected solution


data = scipy.io.loadmat('./dataset/ml-100k/ub_base')  # load 100k dataset

Rating = data['Rating'].flatten()
UserID = data['UserID'].flatten(
) - 1  # Python indexing starts from 0 whereas Matlab from 1
MovID = data['MovID'].flatten(
) - 1  # Python indexing starts from 0 whereas Matlab from 1

nM = np.amax(data['MovID'])
nU = np.amax(data['UserID'])

Z = csr_matrix((Rating, (MovID, UserID)), shape=(nM, nU),
               dtype=float).toarray()
kappa = 5000

tstart = time()
Z_proj = proj_nuc(Z, kappa)
elapsed = time() - tstart
print('proj for 100k data takes {} sec'.format(elapsed))

# NOTE: This one can take few minutes!
data = scipy.io.loadmat('./dataset/ml-1m/ml1m_base')  # load 1M dataset

Rating = data['Rating'].flatten()
UserID = data['UserID'].flatten(
) - 1  # Python indexing starts from 0 whereas Matlab from 1
MovID = data['MovID'].flatten(
) - 1  # Python indexing starts from 0 whereas Matlab from 1

nM = np.amax(data['MovID'])
nU = np.amax(data['UserID'])

Z = csr_matrix((Rating, (MovID, UserID)), shape=(nM, nU),
               dtype=float).toarray()
kappa = 5000

tstart = time()
Z_proj = proj_nuc(Z, kappa)
elapsed = time() - tstart
print('proj for 1M data takes {} sec'.format(elapsed))
