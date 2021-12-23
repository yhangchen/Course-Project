import numpy as np
from scipy.fftpack import fft2, fftshift, ifft2
import matplotlib.pyplot as plt
import scipy.sparse as sps
import cv2
#N.B. This can be deprecated depending on the python version. Use cv2 in this case
#import cv2 # if imread does not work for you.

# N.B. You can install cv2 by calling 'conda install opencv' or 'pip install opencv' in the terminal
# If you have incompatibilities with hdf5, you can uninstall it by command line 'conda remove hdf5' or
# 'pip remove hdf5' and then install opencv

from scipy.sparse.linalg import LinearOperator, svds
from pywt import wavedec2, waverec2, coeffs_to_array, array_to_coeffs
from collections import namedtuple
from math import sqrt

import time
import pylab as pl
from IPython import display

# Mappings between n dimensional complex space and 2n dimensional real space
real2comp = lambda x: x[0:x.shape[0] // 2] + 1j * x[x.shape[0] // 2:]
comp2real = lambda x: np.append(x.real, x.imag)

# Load Data
# x = imread('blurredplate.jpg', flatten=True, mode='F')
x = cv2.imread('blurredplate.jpg',
               cv2.IMREAD_GRAYSCALE)  #If imread does not work for you

x = x[60:188, 40:296]
x = x / np.linalg.norm(x, ord=np.inf)
#x = x/np.linalg.norm(x)
imsize1 = x.shape[0]
imsize2 = x.shape[1]
imsize = imsize1 * imsize2

# Shows the blurred image!
plt.imshow(x, cmap='gray')
plt.show()

# Reshaping operators matrix to vector and vector to matrix
mat = lambda x: np.reshape(x, [imsize1, imsize2])
vec = lambda x: x.flatten()

# Set the measurement vector b
b = comp2real(vec(fft2(fftshift(x))))

# Roughly estimate the support of the blur kernel
K1 = 30
K2 = 30
Indw = np.zeros([imsize1, imsize2])
ind1 = np.int(imsize1 / 2 - (K1 + 1) / 2 + 1)
ind2 = np.int(imsize1 / 2 + (K1 + 1) / 2)
ind3 = np.int(imsize2 / 2 - (K2 + 1) / 2 + 1)
ind4 = np.int(imsize2 / 2 + (K2 + 1) / 2)
Indw[ind1:ind2, ind3:ind4] = 1
#above, for implementational simplicity we assume K1 and K2 odd, even
#if they are even 1 pixel probably won't cause much trouble
plt.imshow(-Indw, cmap='gray')  # Shows the estimated support of blur kernel!
plt.show()
Indw = vec(Indw)
kernelsize = np.count_nonzero(Indw)
Indi = np.nonzero(Indw > 0)[0]
Indv = Indw[Indi]

# Define operators Bop and Cop
Bmat = sps.csr_matrix((Indv, (Indi, range(0, kernelsize))),
                      shape=(imsize, kernelsize))
Bop = lambda x: mat(Bmat.dot(x))
BTop = lambda x: Bmat.T.dot(vec(x))

# Compute and display wavelet coefficients of the original and blurred image
l = coeffs_to_array(wavedec2(x, 'db1', level=4))[1]
Cop = lambda x: waverec2(array_to_coeffs(mat(x), l, output_format='wavedec2'),
                         'db1')
CTop = lambda x: coeffs_to_array(wavedec2(x, 'db1', level=4))[0]

# Define operators
Aoper = lambda m, n, h: comp2real(1.0 / sqrt(imsize) * n * vec(
    fft2(Cop(m)) * fft2(Bop(h))))
AToper = {
    "matvec":
    lambda y, w: CTop(
        np.real(fft2(mat(np.conj(real2comp(y)) * vec(fft2(Bop(w))))))) / sqrt(
            y.shape[0] / 2.0),
    "rmatvec":
    lambda y, w: BTop(np.real(ifft2(mat(real2comp(y) * vec(ifft2(Cop(w))))))) *
    (y.shape[0] / 2.0)**1.5
}


def plot_func(mEst, C, x, iteration):

    xEst = -C(mEst)
    xEst = xEst - min(xEst.flatten())
    xEst = xEst / max(xEst.flatten())
    plt.imshow(xEst, cmap='gray')
    if iteration % 20 == 0:
        plt.savefig(f"figs/{kappa}-{K1}-{K2}-{iteration}.png")
    plt.show()


def frank_wolfe(Aoper, AToper, b, n1, n2, kappa, maxit, plotFunc):
    # PURPOSE: We will solve the following problem formulation with
    # Frank-Wolfe's method.
    #                   min_x  0.5*norm(A(x) - b)^2
    #                   s.t.:  norm_nuc(x) <= kappa,
    #
    # Laboratory for Information and Inference Systems (LIONS)
    # Ecole Polytechnique Federale de Lausanne (EPFL) - SWITZERLAND
    # Last modification: November 26, 2019

    # Print the caption

    #Initialize
    AX_t = 0.0  # zeros
    X = 0.0  # zeros

    # keep track of objective value
    fx = np.array([])

    # The main loop
    for iteration in range(0, maxit + 1):

        # Print the objective values ...
        fx = np.append(fx, 0.5 * np.linalg.norm(AX_t - b, 2)**2)
        print('{:03d} | {:.4e}'.format(iteration, fx[-1]))

        # Form the residual and fix the operator to be used in svds.
        res_cur = AX_t - b
        ATop1 = lambda w: AToper["matvec"](res_cur, w)
        ATop2 = lambda w: AToper["rmatvec"](res_cur, w)
        svdsArg = LinearOperator((n2, n1), matvec=ATop1, rmatvec=ATop2)
        topLe_vec, singVal, topRe_vec = svds(svdsArg,
                                             k=1,
                                             tol=1e-4,
                                             which='LM')

        # Apply A to the rank 1 update
        AXsharp_t = Aoper(topLe_vec, -kappa, topRe_vec.T)

        # Step size
        weight = 2 / (iteration + 2)

        # Update A*X
        AX_t = (1.0 - weight) * AX_t + weight * (AXsharp_t)

        # Update X
        X = (1.0 - weight) * X + weight * (-kappa * topLe_vec @ topRe_vec)

        # Show the reconstruction (at every 10 iteration)
        if (float(iteration / 10) == float(iteration) / 10):
            U, S, V = np.linalg.svd(X, full_matrices=0, compute_uv=1)
            plotFunc(U[:, 0], iteration)

    return X


# Run Frank-Wolfe's method
MaxIters = 200
kappa = 1000
plotF = lambda m, iteration: plot_func(m, Cop, x, iteration)
xFW = frank_wolfe(Aoper, AToper, b, kernelsize, imsize, kappa, MaxIters, plotF)

# NOTE: This experiment is based on the theory and the codes publised in
#'Blind Deconvolution using Convex Programming' by A.Ahmed, B.Recht and J.Romberg.