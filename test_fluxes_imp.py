# This version of test_deblending is much simpler than the original one, because the flux of the weakest source is 
# the same in all ~10000 maps.
import numpy as np
import pyfits as pyf
import pylab as pyl
# import random
from matplotlib.font_manager import fontManager, FontProperties
from scipy import optimize
from scipy import integrate as integr
import scipy.ndimage as nd
import scipy.special.basic as spb
import math
from tkp_lib import utilities, image
import tkp_lib.accessors as access
import tkp_lib.image as imag

outputpath='/home/hspreeuw/validation/fluxes/results/'
n_of_steps=9980
n_of_steps=1000
lowest_flux=10.
highest_flux=1000.
inserted_fluxes=np.linspace(lowest_flux,highest_flux,n_of_steps)
number_of_high_pixels=np.zeros((n_of_steps),'d')
# Keep track of all pixels higher than 6 sigma.

measured_fluxes=np.zeros((n_of_steps,64),'d')
reported_flux_errors=np.zeros((n_of_steps,64),'d')
ra_true=np.zeros((64),'d')
ra_diff=np.zeros((n_of_steps,64),'d')
ra_err=np.zeros((n_of_steps,64),'d')
dec_true=np.zeros((64),'d')
dec_diff=np.zeros((n_of_steps,64),'d')
dec_err=np.zeros((n_of_steps,64),'d')

pathout='/home/hspreeuw/validation/fluxes/FITS/'
edgeskip=16
mapsize=32

degtorad = np.pi / 180.0
raincr_rad = -3.333333414E-03 * degtorad
decincr_rad = 3.333333414E-03 * degtorad
alpha0 = 2.66363244382E+02 * degtorad
delta0 = -2.99529359725E+01 * degtorad
centrerapix=128.
centredecpix=129.

def pix_to_position(xpix, ypix):
        # Converts pixel coordinates to celestial coordinates for a SIN projection map with 
        # tangent point (alpha0,delta0) at pixelcoordinates (centrerapix,centredecpix) and pixelsizes raincr_rad and decincr_rad.
        # The "+ 1.0" in the two lines below is needed because the bottom left pixel has indices [0,0] in numpy.
        # So the central pixel has coordinates [centrerapix-1,centredecpix-1] in numpy.

        L = (xpix + 1.0 - centrerapix) * raincr_rad
        M = (ypix + 1.0 - centredecpix)* decincr_rad

        sqL=np.power(L,2)
        sqM=np.power(M,2)
        help1=np.sqrt(1.0 - sqL - sqM)
        sind=np.sin(delta0)
        cosd=np.cos(delta0)
        help2=M*sind

        alpha=alpha0+np.arctan(L/(help1*cosd-help2))
        delta=np.arcsin(M*cosd+help1*sind)
   
        return alpha/degtorad,delta/degtorad

centrexstart=13.
centreystart=17.
dx=edgeskip
dy=edgeskip
# Now, compute the 64 true ra's and decs.
q=0
# print 'this is the separation index:',1+math.fmod(i-1,number_of_separation_steps)
while q<8:
    p=0
    while p<8:
        xstart=mapsize*q
        # x=xstart+edgeskip
        ystart=mapsize*p
        # y=ystart+edgeskip
        ra_true[8*q+p],dec_true[8*q+p]=pix_to_position(xstart+centrexstart,ystart+centreystart)
        print 'ra_true[8*q+p],dec_true[8*q+p]=',ra_true[8*q+p],dec_true[8*q+p]
        # ra_true[8*q+p]=ra
        # dec_true[8*q+p]=dec
        # sources[x-dx:x+dx,y-dy:y+dy]+=gaussians[0,:,:]*inserted_fluxes[i-1]
        p+=1
    q+=1

# Also, make a grid of intermediate ra,decs to order the measured ra,decs for comparison with the true ra,decs.
# It is handy to make it cover the entire set of values, that is why it has 81 elements instead of 64.
ra_intermediate=np.zeros((81),'d')
dec_intermediate=np.zeros((81),'d')
q=0
while q<9:
    p=0
    while p<9:
        ra_intermediate[9*q+p],dec_intermediate[9*q+p]=pix_to_position(mapsize*q,mapsize*p)
        print 'ra_intermediate[9*q+p],dec_intermediate[9*q+p]=',ra_intermediate[9*q+p],dec_intermediate[9*q+p]
        p+=1
    q+=1
        
def order_positions(ra_array,dec_array):
    # Given 64 Right Ascensions and 64 declinations, actually 64 (ra, dec) combinations,
    # in the form of ra_array((64),'d') and dec_array((64),'d') it returns two
    # ordered arrays of Right Ascensions and declinations.
    ra_ordered=np.zeros((64),'d')
    dec_ordered=np.zeros((64),'d')
    q=0
    while q<8:
        p=0
        while p<8:
            # Right ascension increases from left to right.
            # condition=(ra_array>ra_intermediate[9*q+p+1]) & (ra_array<ra_intermediate[9*q+p])\
            #          & (dec_array<dec_intermediate[9*q+p+1]) & (dec_array>dec_intermediate[9*q+p])
            condition1=(ra_array>ra_intermediate[9*(q+1)+p+1]) & (ra_array<ra_intermediate[9*q+p])
            # print 'ra_intermediate[9*q+p+1],ra_intermediate[9*q+p]=',ra_intermediate[9*q+p+1],ra_intermediate[9*q+p]
            # print 'condition1=',condition1
            condition2=(dec_array<dec_intermediate[9*(q+1)+p+1]) & (dec_array>dec_intermediate[9*q+p])
            # print 'condition2=',condition2
            condition=condition1 & condition2
            # print 'ra_array[condition],dec_array[condition]=',ra_array[condition],dec_array[condition]
            ra_ordered[8*q+p]=ra_array[condition]
            # condition_dec=(dec_array[8*q+p]<dec_intermediate[8*(q+1)+p+1]) & (dec_array[8*q+p]>dec_intermediate[8*q+p])
            dec_ordered[8*q+p]=dec_array[condition]
            p+=1
        q+=1
    return  ra_ordered,dec_ordered
ra_test,dec_test=order_positions(ra_true,dec_true)
print 'This is the result of the reordering of the true ra values,ra_true-ra_test:',ra_true-ra_test
print 'This is the result of the reordering of the true dec values,dec_true-dec_test:',dec_true-dec_test
answer=str(raw_input('Do you want to proceed?'))
ra_true=ra_test
dec_true=dec_test
# Make a series of Gaussians with separations linearly increasing.
# We'll insert those in the maps. They are separated perpendicular to the major axis, so along the minor axis.
number_of_separation_steps=20
gaussians=np.zeros((number_of_separation_steps+1,mapsize,mapsize),'d')
centrex=np.zeros((number_of_separation_steps+1),'d')
centrey=np.zeros((number_of_separation_steps+1),'d')
pixelsize=12.
theta=-np.pi*49.80/180.
semimajor=0.5*67.1472/pixelsize
semiminor=0.5*56.1528/pixelsize
# Try with circularly symmetric Gaussian
# semimajor=semiminor

minor_width=semiminor/math.sqrt(2.*math.log(2.))
print 'minor_width=',minor_width
beginsep=2.0*minor_width
endsep=7.0*minor_width

j=0
while j<=number_of_separation_steps:
    # I conveniently use sign(j) because sign(0)=0 and sign(1,2,3..)=1.
    centrex[j]=centrexstart+np.sign(j)*np.cos(theta)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
    centrey[j]=centreystart+np.sign(j)*np.sin(theta)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
    print 'this is j and the separation:',j,np.sign(j)*(beginsep+float(j-1)*(endsep-beginsep)/float(number_of_separation_steps-1))
    k=0
    while k<mapsize:
        x=float(k)
        l=0
        while l<mapsize:
            y=float(l)
            gaussians[j,k,l]=np.exp(-math.log(2.)*(((np.cos(theta)*(x-centrex[j])+np.sin(theta)*(y-centrey[j]))/semiminor)**2.+((np.cos(theta)*(y-centrey[j])-np.sin(theta)*(x-centrex[j]))/semimajor)**2.))
            l+=1
        k+=1
    j+=1


central_pos=nd.center_of_mass(gaussians[0,:,:])
print 'central_pos=',central_pos[0],central_pos[1]
# Start loop
number_too_much_noise=0
number_not_64=0
i=1
while i<=n_of_steps:
 # Import the noise map.
 hdulist=pyf.open('/home/hspreeuw/validation/FITS2/FITS/PURENOISE'+str(i)+'.FITS')
 dummyhigh=np.where(hdulist['PRIMARY'].data>lowest_flux,1,0)
 number_of_high_pixels[i-1]=np.sum(dummyhigh)
 if number_of_high_pixels[i-1]>0:
     print 'there is at least one noise pixel higher than the lowest flux!!!!!!!!!!!!!!!!!!!'
     number_too_much_noise+=1
 print 'number too much noise=',number_too_much_noise
 print 'i=',i 
 
 sources=np.zeros((256,256),'d')

 q=0
 # print 'this is the separation index:',1+math.fmod(i-1,number_of_separation_steps)
 while q<8:
     p=0
     while p<8:
         xstart=mapsize*q
         x=xstart+edgeskip
         ystart=mapsize*p
         y=ystart+edgeskip
         sources[x-dx:x+dx,y-dy:y+dy]+=gaussians[0,:,:]*inserted_fluxes[i-1]
         p+=1
     q+=1
 hdulist['PRIMARY'].data+=sources.transpose()
 scidata=hdulist['PRIMARY'].data
 # These should be exact copies of the source free maps!!!
 hdulist.writeto(pathout+'SOURCESINSERTED'+str(i)+'.FITS')
 hdulist.close()
 
 # The file is opened immediately after it has been made.
 my_fitsfile = access.FitsFile(pathout+'SOURCESINSERTED'+str(i)+'.FITS')
 my_image= imag.ImageData(my_fitsfile)
 # For the source extraction, both the detection and the analysis threshold are set equal to lowerlim in settings.py.
 sextract_results = my_image.sextract()

 if len(sextract_results)!=64:
     print 'len(sextract_results)=',len(sextract_results)
 else:
     # Only leap to the next image if there were exactly 64 detections.
     ra_values=np.zeros((64),'d')
     dec_values=np.zeros((64),'d')
     k=0
     for detection in sextract_results:
         measured_fluxes[i-1,k]=detection.peak.value
         reported_flux_errors[i-1,k]=detection.peak.error
         ra_values[k]=detection.ra.value
         dec_values[k]=detection.dec.value
         ra_err[i-1,k]=detection.ra.error
         dec_err[i-1,k]=detection.dec.error
         # print 'detection.ra,detection.dec=',detection.ra,detection.dec
         k+=1
     ra_values,dec_values=order_positions(ra_values,dec_values)
     # print 'ra_values=',ra_values
     # print 'ra_true=',ra_true
     ra_diff[i-1,:]=ra_values-ra_true
     dec_diff[i-1,:]=dec_values-dec_true
     print 'np.std(ra_diff[i-1,:]),np.std(dec_diff[i-1,:])=',np.std(ra_diff[i-1,:]),np.std(dec_diff[i-1,:])
 number_not_64+=len(sextract_results)
 print 'number_not_64,np.mod(number_not_64,64)=',number_not_64,np.mod(number_not_64,64)
 i+=1

np.save(outputpath+'inserted_fluxes',inserted_fluxes)
np.save(outputpath+'measured_fluxes',measured_fluxes)
np.save(outputpath+'reported_flux_errors',reported_flux_errors)
np.save(outputpath+'number_of_high_pixels',number_of_high_pixels)
np.save(outputpath+'ra_diff',ra_diff)
np.save(outputpath+'ra_err',ra_err)
np.save(outputpath+'dec_diff',dec_diff)
np.save(outputpath+'dec_err',dec_err)

