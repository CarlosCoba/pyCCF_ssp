import numpy as np
speed_light=3e5 #km/s

def set_parameters(z=None):
	#elines to fit

	elines=[6548.1001,6562.8169,6583.6001]
	#elines=[4958.9102,5006.8398,np.nan]

	if z==None:
	#redshift
		z=0.000
		redshift=[z,z-400./speed_light,z+400./speed_light]
	else:
		redshift=[z,z-400./speed_light,z+400./speed_light]		


	#FIX PARAMETER? 0:NO 1:YES,value,min,max
	sigma=[0,1.1,1,np.inf]

	#FIT AMPLITUDE?	0:NO 1:YES
	link=[0,1,1]

	#ELINE TO LINK
	eline_to_link=[2,None,None]

	# MULTIPLICATIVE LINK
	type=[1/3.,1,1]

	return [elines,redshift,sigma,link,eline_to_link,type]


