#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np
import pyfits
from astropy.io import fits
from PyAstronomy import pyasl
from lmfit import Minimizer, Parameters, report_fit
import sys
#import scipy.interpolate as sci
from scipy.signal import square, sawtooth, correlate
import scipy
from r_ccf import r_ccf
import scipy.interpolate as sci






nargs=len(sys.argv)

def main_program():
	if (nargs==7):
		cube=sys.argv[1]
		template_ssp=sys.argv[2]
		segmentation_map=sys.argv[3]
		z=float(sys.argv[4])
		name=sys.argv[5]
		plot=int(sys.argv[6])
		return cube,template_ssp,segmentation_map,z,name,plot
		
	else:
		print 'USE: CCF_IFS.py cube.fits template_ssp segmentation_map z OUTPUT plot'
		exit()




speed_light=299792.485#km/s


def ccf(cube,template_ssp,segmentation_map,z,name,plot):
	#import this scripts

	data, head = fits.getdata(cube, header=True,memmap=True)
	wave = head['CRVAL3'] + head['CDELT3'] * np.arange(head['NAXIS3'])
	#wave = head['CRVAL3'] + head['CD3_3'] * np.arange(head['NAXIS3'])
	print np.min(wave)


	[nz,ny,nx]=data.shape
######################
	# template ssp
	#rss=0
	hdu_temp = fits.open(template_ssp)
	data_temp = hdu_temp[0].data
	head_temp = hdu_temp[0].header             
	wave_tmp = head_temp['CRVAL1'] + head_temp['CDELT1'] * np.arange(head_temp['NAXIS1']) 
	print np.min(wave_tmp)
	#org_spec=data_temp[0,rss,:]
	#model_spec=data_temp[1,rss,:]
	#model_joint_spec=data_temp[2,rss,:]
	[nz_tmp,ny_tmp,nx_temp]=data_temp.shape
######################
	# cont segmentation
	data_seg=fits.getdata(segmentation_map)


	#Remove NANs from cube:
	data[np.isnan(data)]=0

	#Prepare arrays to save data
	MEDIAN=np.ones((ny,nx))*np.nan
	MEAN=np.ones((ny,nx))*np.nan
	BISECTOR=np.ones((9,ny,nx))*np.nan
	RAD_VEL=np.ones((ny,nx))*np.nan
	DISPERSION=np.ones((ny,nx))*np.nan
	DISPERSION_ELINE=np.ones((ny,nx))*np.nan

	for jj in range(ny):
		for ii in range(nx):
			try:
					print jj,ii


					# the bin ==0 in the segmentation map is not a real bin:
					seg_val=int(data_seg[jj][ii])-1
					seg_val=0
					if seg_val<0:
						seg_val=np.nan
					model_spec=data_temp[1,seg_val,:]
					

					Mgb=5172*(1+z)
					#NaD=5892*(1+z)
					#Mgb=NaD
					index=np.where((wave>Mgb-50) & (wave<Mgb+50))
					index2=np.where((wave>Mgb-100) & (wave<Mgb+100))
				
					flux=data[:,jj,ii]
					flux=flux[index]
					flux_obs=flux/np.max(flux)
					wave_obs=wave[index]
					wave_temp=wave_tmp[index]


					flux_temp=model_spec[index]
					flux_temp=flux_temp/np.nanmax(flux_temp)







					###########################################################
					#This is the real cross correlation function with 20 lags:#
					###########################################################
					lags_true,r_true=r_ccf(flux_temp,flux,True,20)
					




					flux_temp2=model_spec[index2]
					flux_temp2=flux_temp2/np.max(flux_temp2)



					rv=np.arange(-700,700,5)
					ccf_vel=[]
					for j,delta_v in enumerate(rv):
						fi = sci.interp1d(wave_tmp[index2]*(1.0 + delta_v/speed_light), flux_temp2)
						# The values to interpolate must be in the range covered by xt*(1.0 + delta_v/speed_light)
						lags,r_=r_ccf(fi(wave_obs),flux_obs,False,1)
						ccf_vel.append(r_[0])



					# The cross correlation in velocity space:
					ccf_vel=np.asarray(ccf_vel)
					# Normalize the ccf:
					ccf_vel=ccf_vel/np.nanmax(ccf_vel)

					# Estimate the peak
					max_ccf_index=ccf_vel.argmax()
					max_ccf=ccf_vel[max_ccf_index]
					mean_vel=rv[max_ccf_index]
					v_peak=mean_vel



					# PERFORM A QUICK FIT TO ESTIMATE VELOCITY AND SIGMA 
					# AROUND THE PEAK
					from fit_gaussians import G1
					best_fit_ccf=G1(rv[rv>v_peak-500],ccf_vel[rv>v_peak-500],max_ccf,v_peak,100,0)
					A_ccf,v_ccf,sigma_ccf=best_fit_ccf[0],best_fit_ccf[1],abs(best_fit_ccf[2])

					v_centroid=v_ccf
					#shift velocity to the centroid
					rv=rv-v_centroid


					# shift cross correlation to the region of inrerenst -450km/s<rv<450 km/s
					index_vel=np.where((rv>-450) & (rv<450))
					rv=rv[index_vel]
					ccf_vel=ccf_vel[index_vel]

					#Calculate the bisectors
					x_bis=[]
					y_bis=[]
					delta_v=[]
					step=0.1
					for i in np.arange(step,1,step):
						li=(1-i)
						ls=(1+step-i)
						index=np.where((ccf_vel>li) & (ccf_vel<=ls))
						vels=np.mean(rv[index])
						x_bis.append(vels)
						y_bis.append(li)
						delta_v.append(vels)


					#SAVE CCF and BISECTORS
					J,I= int(jj),int(ii)
					MEAN[J][I]= np.nanmean(delta_v)
					MEDIAN[J][I]=np.nanmedian(delta_v)
					BISECTOR[:,J,I]=delta_v





					#SAVE VELOCITY PARAMETERS FROM FIT
					J,I= int(jj),int(ii)
					RAD_VEL[J][I]= v_ccf+0
					DISPERSION[J][I]=sigma_ccf+0
					DISPERSION_ELINE[J][I]=0


					if plot == 1:	
							plt.figure(figsize=(8,4))

							plt.subplot(131)
							plt.plot(wave_temp,flux_temp,label="template")
							plt.plot(wave_obs,flux_obs,'k-',label="spectra")
							plt.xlabel('wavelength')
							plt.legend(loc="upper left")

							plt.subplot(132)
							plt.plot(lags_true,r_true,"r-",label='r_ccf')
							plt.axhline(y=0, color='r', linestyle='--')
							plt.xlabel('wavelength')
							plt.ylabel('r_ccf')
							plt.ylim(-1,1)
							plt.legend(loc="upper left")

							plt.subplot(133)
							plt.plot(rv,ccf_vel,'k+',label='norm. ccf_vel')
							G1(rv,ccf_vel,max_ccf,mean_vel,100,1)
							plt.xlabel('velocity')
							plt.ylabel('r_ccf')
							plt.axhline(y=0, color='r', linestyle='--')
							plt.ylim(-1,1.5)
							plt.plot(x_bis,y_bis,'g.')
							plt.legend(loc="upper left")
							plt.tight_layout()

							plt.show()


			except (ValueError,ZeroDivisionError,IOError,TypeError,IndexError):
			#except (1):
				pass




	###
	# save all
	###

	hdu_1 = fits.PrimaryHDU()
	hdu_1.data=MEDIAN
	hdu_1.header=head
	hdu_1.writeto("%s.ccf.median.fits"%name)


	hdu_2 = fits.PrimaryHDU()
	hdu_2.data=MEAN
	hdu_2.header=head
	hdu_2.writeto("%s.ccf.mean.fits"%name)


	hdu_3 = fits.PrimaryHDU()
	hdu_3.data=BISECTOR
	hdu_3.header=head
	hdu_3.writeto("%s.ccf.bisector.cube.fits"%name)

	hdu_4 = fits.PrimaryHDU()
	hdu_4.data=RAD_VEL
	hdu_4.header=head
	hdu_4.writeto("%s.ccf.vel.fits"%name)

	hdu_5 = fits.PrimaryHDU()
	hdu_5.data=DISPERSION
	hdu_5.header=head
	hdu_5.writeto("%s.ccf.sigma_ccf.fits"%name)

	hdu_6 = fits.PrimaryHDU()
	hdu_6.data=DISPERSION_ELINE
	hdu_6.header=head
	hdu_6.writeto("%s.ccf.sigma_eline.fits"%name)


	print "done"

if __name__ == "__main__":
	main=main_program()
	ccf(main[0],main[1],main[2],main[3],main[4],main[5])
#else:
#	ccf()	







