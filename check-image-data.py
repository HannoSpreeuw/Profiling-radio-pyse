from astropy.io import fits
import numpy as np

input_fits = 'SOURCESINSERTED_100Jy.FITS'
hdulist=fits.open(input_fits)
imagedata=hdulist['PRIMARY'].data
dimensions = (hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2'])

average=np.mean(imagedata)
min=np.min(imagedata)
max=np.max(imagedata)
std=np.std(imagedata)

print()
print("The dimensions are {}".format(dimensions))
print("mean={}, min={}, max={}, std={}".format(average, min, max, std))

random_data=np.float32(np.random.normal(size=dimensions))
average=np.mean(random_data)
min=np.min(random_data)
max=np.max(random_data)
std=np.std(random_data)
print()
print("Random data of these dimensions have these properties:")
print("mean={}, min={}, max={}, std={}".format(average, min, max, std))
