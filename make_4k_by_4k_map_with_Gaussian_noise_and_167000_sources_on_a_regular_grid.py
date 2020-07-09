# This version of test_deblending is much simpler than the original one, because the flux of the weakest source is 
# the same in all ~10000 maps.
import numpy as np
from astropy.io import fits
import pylab as pyl
import os

input_fits = '89M18BIGPBC.FITS'
output_fits = 'SOURCESINSERTED.FITS'

degtorad = np.pi / 180.0
raincr_rad = -3.333333414E-03 * degtorad
decincr_rad = 3.333333414E-03 * degtorad
alpha0 = 2.66363244382E+02 * degtorad
delta0 = -2.99529359725E+01 * degtorad


# # Make a series of Gaussians with separations linearly increasing.
# # We'll insert those in the maps. They are separated perpendicular to the major axis, so along the minor axis.
# number_of_separation_steps=20
# gaussians=np.zeros((number_of_separation_steps+1,mapsize,mapsize),'d')
# centrex=np.zeros((number_of_separation_steps+1),'d')
# centrey=np.zeros((number_of_separation_steps+1),'d')
# pixelsize=12.
# theta=-np.pi*49.80/180.
# semimajor=0.5*67.1472/pixelsize
# semiminor=0.5*56.1528/pixelsize
# 
# minor_width=semiminor/math.sqrt(2.*math.log(2.))
# print 'minor_width=',minor_width
# beginsep=2.0*minor_width
# endsep=7.0*minor_width
# 
# j=0
# while j<=number_of_separation_steps:
#     # I conveniently use sign(j) because sign(0)=0 and sign(1,2,3..)=1.
#     centrex[j]=centrexstart+np.sign(j)*np.cos(theta)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
#     centrey[j]=centreystart+np.sign(j)*np.sin(theta)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
#     print 'this is j and the separation:',j,np.sign(j)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
#     k=0
#     while k<mapsize:
#         x=float(k)
#         l=0
#         while l<mapsize:
#             y=float(l)
#             gaussians[j,k,l]=np.exp(-math.log(2.)*(((np.cos(theta)*(x-centrex[j])+np.sin(theta)*(y-centrey[j]))/semiminor)**2.+((np.cos(theta)*(y-centrey[j])-np.sin(theta)*(x-centrex[j]))/semimajor)**2.))
#             l+=1
#         k+=1
#     j+=1


hdulist=fits.open(input_fits)
dimensions = (hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2'])
print
print "The dimensions of this image are: ",dimensions

sources=np.zeros(dimensions,'d')

# q=0
# # print 'this is the separation index:',1+math.fmod(i-1,number_of_separation_steps)
# while q<8:
#     p=0
#     while p<8:
#         xstart=mapsize*q
#         x=xstart+edgeskip
#         ystart=mapsize*p
#         y=ystart+edgeskip
#         sources[x-dx:x+dx,y-dy:y+dy]+=gaussians[0,:,:]*inserted_fluxes[i-1]
#         p+=1
#     q+=1

# This gives us random noise, normally distributed, as single precision floats.
hdulist['PRIMARY'].data = np.float32(np.random.normal(size=dimensions))
scidata=hdulist['PRIMARY'].data
# These should be exact copies of the source free maps!!!
try: 
    hdulist.writeto(output_fits)
except IOError:
    print "Output file already exists, replacing"
    os.remove(output_fits)
    hdulist.writeto(output_fits)

hdulist.close()
 
#  # The file is opened immediately after it has been made.
#  my_fitsfile = access.FitsFile(pathout+'SOURCESINSERTED'+str(i)+'.FITS')
#  my_image= imag.ImageData(my_fitsfile)
#  # For the source extraction, both the detection and the analysis threshold are set equal to lowerlim in settings.py.
#  sextract_results = my_image.sextract()
# 
#  if len(sextract_results)!=64:
#      print 'len(sextract_results)=',len(sextract_results)
#  else:
#      # Only leap to the next image if there were exactly 64 detections.
#      ra_values=np.zeros((64),'d')
#      dec_values=np.zeros((64),'d')
#      k=0
#      for detection in sextract_results:
#          measured_fluxes[i-1,k]=detection.peak.value
#          reported_flux_errors[i-1,k]=detection.peak.error
#          ra_values[k]=detection.ra.value
#          dec_values[k]=detection.dec.value
#          ra_err[i-1,k]=detection.ra.error
#          dec_err[i-1,k]=detection.dec.error
#          # print 'detection.ra,detection.dec=',detection.ra,detection.dec
#          k+=1
#      ra_values,dec_values=order_positions(ra_values,dec_values)
#      # print 'ra_values=',ra_values
#      # print 'ra_true=',ra_true
#      ra_diff[i-1,:]=ra_values-ra_true
#      dec_diff[i-1,:]=dec_values-dec_true
#      print 'np.std(ra_diff[i-1,:]),np.std(dec_diff[i-1,:])=',np.std(ra_diff[i-1,:]),np.std(dec_diff[i-1,:])
#  number_not_64+=len(sextract_results)
#  print 'number_not_64,np.mod(number_not_64,64)=',number_not_64,np.mod(number_not_64,64)
#  i+=1
# 
# np.save(outputpath+'inserted_fluxes',inserted_fluxes)
# np.save(outputpath+'measured_fluxes',measured_fluxes)
# np.save(outputpath+'reported_flux_errors',reported_flux_errors)
# np.save(outputpath+'number_of_high_pixels',number_of_high_pixels)
# np.save(outputpath+'ra_diff',ra_diff)
# np.save(outputpath+'ra_err',ra_err)
# np.save(outputpath+'dec_diff',dec_diff)
# np.save(outputpath+'dec_err',dec_err)

