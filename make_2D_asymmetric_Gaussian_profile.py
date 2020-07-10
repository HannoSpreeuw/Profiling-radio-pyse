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

def gkern_asym(kernlen=10, peak=10, stdx=1, stdy=1.5):
    """Returns a 2D Gaussian kernel array."""
    gkernx = signal.gaussian(kernlen, std=stdx)
    gkerny = signal.gaussian(kernlen, std=stdy)
    gkern2d = np.outer(gkernx, gkerny)
    return peak*gkern2d

asym_gauss = gkern_asym()
pyl.imshow(asym_gauss)
pyl.show()

all_sources=np.tile(asym_gauss,(409,409))
pyl.imshow(all_sources)
pyl.show()
