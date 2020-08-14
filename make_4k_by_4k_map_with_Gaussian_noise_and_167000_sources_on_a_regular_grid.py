# This version of test_deblending is much simpler than the original one, because the flux of the weakest source is 
# the same in all ~10000 maps.
import numpy as np
import scipy.stats as st
from astropy.io import fits
from scipy import signal
import pylab as pyl
import os, math

input_fits = '89M18BIGPBC.FITS'
output_fits = 'SOURCESINSERTED.FITS'

number_of_sources = 167000

hdulist=fits.open(input_fits)
dimensions = (hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2'])
print()
print("The dimensions of this image are: ",dimensions)

all_sources_expanded_to_image_size=np.zeros(dimensions,'d')

number_of_sources_along_equal_axes = int(math.ceil(np.sqrt(number_of_sources)))
number_of_sources_along_xaxis = int(number_of_sources_along_equal_axes*dimensions[0]/dimensions[1])
number_of_sources_along_yaxis = int(number_of_sources_along_equal_axes*dimensions[1]/dimensions[0])

print()
print("number_of_sources_along_xaxis = ",number_of_sources_along_xaxis,"number_of_sources_along_yaxis = ",number_of_sources_along_yaxis)

kernlenx = int(dimensions[0]/number_of_sources_along_xaxis)
kernleny = int(dimensions[1]/number_of_sources_along_yaxis)
print()
print("kernlenx = ",kernlenx,"  kernleny = ",kernleny)

# Here I try to make enough copies of a 2D asymmetrical Gaussian to cover the entire map, with peaks well above the noise
def gkern_asym(kernlenx, kernleny, peak=100, stdx=1, stdy=1.5):
    """Returns a 2D Gaussian kernel array."""
    gkernx = signal.gaussian(kernlenx, std=stdx)
    gkerny = signal.gaussian(kernleny, std=stdy)
    gkern2d = np.outer(gkernx, gkerny)
    return peak*gkern2d

asym_gauss = gkern_asym(kernlenx, kernleny)
all_sources=np.tile(asym_gauss,(number_of_sources_along_xaxis, number_of_sources_along_yaxis))

edgex = int((dimensions[0]-all_sources.shape[0])/2)
edgey = int((dimensions[1]-all_sources.shape[1])/2)
print()
print("edgex = ",edgex,"  edgey = ",edgey)

all_sources_expanded_to_image_size[edgex:edgex+all_sources.shape[0], edgey:edgey+all_sources.shape[1]] = all_sources

# This gives us random noise, normally distributed, as single precision floats with strong sources on a regular grid added to that.
hdulist['PRIMARY'].data = np.float32(np.random.normal(size=dimensions))+all_sources_expanded_to_image_size
# These should be exact copies of the source free maps!!!
try: 
    hdulist.writeto(output_fits)
except IOError:
    print("Output file already exists, replacing")
    os.remove(output_fits)
    hdulist.writeto(output_fits)

hdulist.close()
 
