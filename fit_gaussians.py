from lmfit import Minimizer, Parameters, report_fit
import numpy as np
import matplotlib.pylab as plt
from set_params import set_parameters
from astropy.modeling.models import custom_model

from astropy.stats import sigma_clip
speed_light=3e5 #km/s

## GAUSSIANS TO MODEL EMISSION LINES:

def G1(wave_range,flux_range,A0,eline0,sigma0,plot=0):
			def G_model(x,A0=1,eline0=1,sigma0=300):

				l_0=eline0
				model0=A0*np.exp(-0.5*(x-l_0)**2/(sigma0**2))
				
				return model0


			def G_deriv(x,A0=1,eline0=1,sigma0=300):
				# Jacobian of G_model
				l_0=eline0
				y0=(x-l_0)/(sigma0)
				model0=A0*np.exp(-0.5*y0**2)


				d_A0=np.exp(-0.5*y0**2)
				d_sigma0=A0*d_A0*(x-l_0)**2/sigma0**3
 				d_l0=A0*d_A0*(x-l_0)/sigma0**2

				return [d_A0,d_l0,d_sigma0]




			# initialize fitters
			from astropy.modeling import models, fitting
			fit = fitting.LevMarLSQFitter()
			or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=1, sigma=3.0)
			GaussModel = custom_model(G_model, fit_deriv=G_deriv)
			model=GaussModel()
			or_fitted_model, mask = or_fit(model, wave_range,flux_range)
			#or_fitted_model = fit(model, wave_range,flux_range)

			A0_best,sigma_best,eline_best=mask.A0,mask.sigma0,mask.eline0



			if plot ==1:

				plt.plot(wave_range, flux_range, 'k+')
				plt.plot(wave_range,G_model(wave_range,A0_best,eline_best,sigma_best),"b-",label='best_fit')
				plt.legend()



			return [A0_best,eline_best,sigma_best,or_fitted_model]#best_vals




def fitting(wave_range,flux_range,z,plot=0):
	set_params=set_parameters(z)
	elines=set_params[0]
	redshift=set_params[1]
	z=redshift[0]
	z_min,z_max=redshift[1],redshift[2]
	disp=set_params[2]
	fit_sigma,sigma,sigma_min,sigma_max=disp[0],disp[1],disp[2],disp[3]
	link=np.asarray(set_params[3])
	eline_link=set_params[4]
	type=set_params[5]

	n_elines=np.count_nonzero(~np.isnan(elines))


	if n_elines==1:
			def G_model(x,A0=1,sigma0=2,z00=z):

				l_0=elines[0]*(1+z00)
				model0=A0*np.exp(-0.5*(x-l_0)**2/(sigma0**2))
				
				return model0+model1


			def G_deriv(x,A0=1,sigma0=2,z00=z):
				# Jacobian of G_model
				l_0=elines[0]*(1+z00)

				y0=(x-l_0)/(sigma0)


				model0=A0*np.exp(-0.5*y0**2)


				d_A0=np.exp(-0.5*y0**2)
				d_sigma0=A0*d_A0*(x-l_0)**2/sigma0**3
				d_z00=elines[0]*A0*d_A0*(x-l_0)/sigma0**2

				return [d_A0,d_sigma0,d_z00]



			# initialize fitters
			from astropy.modeling import models, fitting
			fit = fitting.LevMarLSQFitter()
			or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=3, sigma=3.0)
			GaussModel = custom_model(G_model, fit_deriv=G_deriv)
			model=GaussModel()
			or_fitted_model, mask = or_fit(model, wave_range,flux_range)
			A0_best,sigma_best,z_best=mask.A0,mask.sigma0,mask.z00

			if plot ==1:

				plt.plot(wave_range, flux_range, 'k+')
				plt.plot(wave_range,G_model(wave_range,A0_best,sigma_best,z_best),"k-")



			return [sigma_best,[A0_best],z_best,or_fitted_model]#best_vals




	if n_elines==2:
			def G_model(x,A0=1,A1=1,sigma0=2,z00=z):
				l_0=elines[0]*(1+z00)
				l_1=elines[1]*(1+z00)


				model0=A0*np.exp(-0.5*(x-l_0)**2/(sigma0**2))
				model1=A1*np.exp(-0.5*(x-l_1)**2/(sigma0**2))
				
				return model0+model1


			def G_deriv(x,A0=1,A1=1,sigma0=2,z00=z):
				# Jacobian of G_model
				l_0=elines[0]*(1+z00)
				l_1=elines[1]*(1+z00)

				y0=(x-l_0)/(sigma0)
				y1=(x-l_1)/(sigma0)

				model0=A0*np.exp(-0.5*y0**2)
				model1=A1*np.exp(-0.5*y1**2)

				d_A0=np.exp(-0.5*y0**2)
				d_A1=np.exp(-0.5*y1**2)
				d_sigma0=A0*d_A0*(x-l_0)**2/sigma0**3+A1*d_A1*(x-l_1)**2/sigma0**3
				d_z00=elines[0]*A0*d_A0*(x-l_0)/sigma0**2+elines[1]*A1*d_A1*(x-l_1)/sigma0**2

				return [d_A0,d_A1,d_sigma0,d_z00]


			if 0 in link:

				#Amplitud que no se va a ajustar
				link_index_0=np.where(link==0)[0][0]
				#Las que si se van a ajustar
				link_index_1=np.where(link==1)[0][0]
				#A que linea se va a linkear.
				eline_to_link=eline_link[link_index_0]



				if link_index_0==0 and eline_to_link ==1:
					A0_best=A1_best/3.
				if link_index_0==1 and eline_to_link ==0:
					A1_best=A0_best/3.



			# initialize fitters
			from astropy.modeling import models, fitting
			fit = fitting.LevMarLSQFitter()
			or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=3, sigma=3.0)
			GaussModel = custom_model(G_model, fit_deriv=G_deriv)
			model=GaussModel()
			or_fitted_model, mask = or_fit(model, wave_range,flux_range)
			A0_best,A1_best,sigma_best,z_best=mask.A0,mask.A1,mask.sigma0,mask.z00

			if plot ==1:

				plt.plot(wave_range, flux_range, 'k+')
				plt.plot(wave_range,G_model(wave_range,A0_best,A1_best,sigma_best,z_best),"k-")



			return [sigma_best,[A0_best,A1_best],z_best,or_fitted_model]#best_vals

	if n_elines==3:

		################################################################################






			def G_model(x,A0=1,A1=1,A2=1,sigma0=2,z00=z):
				l_0=elines[0]*(1+z00)
				l_1=elines[1]*(1+z00)
				l_2=elines[2]*(1+z00)

				model0=A0*np.exp(-0.5*(x-l_0)**2/(sigma0**2))
				model1=A1*np.exp(-0.5*(x-l_1)**2/(sigma0**2))
				model2=A2*np.exp(-0.5*(x-l_2)**2/(sigma0**2))
				#A0=A2/3.
				return model0+model1+model2


			def G_deriv(x,A0=1,A1=1,A2=1,sigma0=2,z00=z):
				# Jacobian of G_model
				l_0=elines[0]*(1+z00)
				l_1=elines[1]*(1+z00)
				l_2=elines[2]*(1+z00)
				y0=(x-l_0)/(sigma0)
				y1=(x-l_1)/(sigma0)
				y2=(x-l_2)/(sigma0)
				#A0=A2/3.

				model0=A0*np.exp(-0.5*y0**2)
				model1=A1*np.exp(-0.5*y1**2)
				model2=A2*np.exp(-0.5*y2**2)

				d_A0=np.exp(-0.5*y0**2)
				d_A1=np.exp(-0.5*y1**2)
				d_A2=np.exp(-0.5*y2**2)
				d_sigma0=A0*d_A0*(x-l_0)**2/sigma0**3+A1*d_A1*(x-l_1)**2/sigma0**3+A2*d_A2*(x-l_2)**2/sigma0**3
				d_z00=elines[0]*A0*d_A0*(x-l_0)/sigma0**2+elines[1]*A1*d_A1*(x-l_1)/sigma0**2+elines[2]*A2*d_A2*(x-l_2)/sigma0**2


				return [d_A0,d_A1,d_A2,d_sigma0,d_z00]




			# initialize fitters
			from astropy.modeling import models, fitting
			fit = fitting.LevMarLSQFitter()
			or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,niter=1, sigma=3)
			GaussModel = custom_model(G_model, fit_deriv=G_deriv)
			model=GaussModel(bounds={"sigma0":(1.,10.),"z00":(z-400./3e5,z+400./3e5),"A1":(0.2,1.5)})
			or_fitted_model, mask = or_fit(model, wave_range,flux_range)
			fitted_model = fit(model, wave_range,flux_range)
			A0_best,A1_best,A2_best,sigma_best,z_best=mask.A2,mask.A1,mask.A2,mask.sigma0,mask.z00
			A0_best,A1_best,A2_best,sigma_best,z_best=fitted_model.A2.value,fitted_model.A1.value,fitted_model.A2.value,fitted_model.sigma0.value,fitted_model.z00.value










			if 0 in link:

				#Amplitud que no se va a ajustar
				link_index_0=np.where(link==0)[0][0]
				#Las que si se van a ajustar
				link_index_1=np.where(link==1)[0][0]
				#A que linea se va a linkear.
				eline_to_link=eline_link[link_index_0]

				if link_index_0==0 and eline_to_link ==2:
					A0_best= A2_best/3.
				if link_index_0==0 and eline_to_link ==1:
					A0_best= A1_best/3.
				if link_index_0==1 and eline_to_link ==0:
					A1_best= A0_best/3.
				if link_index_0==1 and eline_to_link ==2:
					A1_best= A2_best/3.
				if link_index_0==2 and eline_to_link ==0:
					A2_best= A0_best/3.
				if link_index_0==2 and eline_to_link ==1:
					A2_best= A1_best/3.


		################################################################################
	
			if plot ==1:

				plt.plot(wave_range, flux_range, 'k+')
				plt.plot(wave_range,G_model(wave_range,A0_best,A1_best,A2_best,sigma_best,z_best),"k-")
				plt.plot(wave_range,or_fitted_model,"ro")



			return [sigma_best,[A0_best,A1_best,A2_best],z_best,or_fitted_model]#best_vals




