import numpy as np
import scipy.stats as st
import pylab as pyl

def gkern_sym(kernlen=21, nsig=3):
    """Returns a 2D Gaussian kernel."""

    x = np.linspace(-nsig, nsig, kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kern2d = np.outer(kern1d, kern1d)
    return kern2d/kern2d.sum()

pyl.imshow(gkern_sym())
pyl.show()

from scipy import signal

def gkern_asym(kernlen=21, stdx=2, stdy=3):
    """Returns a 2D Gaussian kernel array."""
    gkernx = signal.gaussian(kernlen, std=stdx)
    gkerny = signal.gaussian(kernlen, std=stdy)
    gkern2d = np.outer(gkernx, gkerny)
    return gkern2d

pyl.imshow(gkern_asym())
pyl.show()
