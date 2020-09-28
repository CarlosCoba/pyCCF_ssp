import numpy as np
import matplotlib.pylab as plt
import scipy.interpolate as sci


def r_ccf(template,data,both=True,nlags=20):

	data=data/np.nanmax(data)
	template=template/np.nanmax(template)

	lag_arr=[]
	r_arr=[]
	n=len(data)

	if both== True:
		for lag in range(-nlags,nlags):
			m=float(n)-abs(lag)
			if lag>=0:

				x_template,x_data=template[:n-lag],data[lag:]
				sum_XY=np.sum(x_template*x_data)
				sum_X,sum_Y=np.sum(x_template),np.sum(x_data)
				cov=(m*sum_XY-sum_X*sum_Y)/(m*(m-1))

				X_square=np.sum(x_template**2)
				Y_square=np.sum(x_data**2)
				deno1=m*Y_square-np.sum(x_data)**2
				deno2=(m*X_square-np.sum(x_template)**2)
				vary=np.sqrt(deno2*deno1)/(m*(m-1))
				r_lag=cov/vary
				lag_arr.append(lag)
				r_arr.append(r_lag)


			if lag<0:
				lag=abs(lag)
				x_template,x_data=data[:n-lag],template[lag:]
				sum_XY=np.sum(x_template*x_data)
				sum_X,sum_Y=np.sum(x_template),np.sum(x_data)
				cov=(m*sum_XY-sum_X*sum_Y)/(m*(m-1))

				X_square=np.sum(x_template**2)
				Y_square=np.sum(x_data**2)
				deno1=m*Y_square-np.sum(x_data)**2
				deno2=(m*X_square-np.sum(x_template)**2)
				vary=np.sqrt(deno2*deno1)/(m*(m-1))
				r_lag=cov/vary
				lag_arr.append(-lag)
				r_arr.append(r_lag)

	else:
		for lag in range(nlags):
			m=float(n)-abs(lag)
			if lag>=0:

				x_template,x_data=template[:n-lag],data[lag:]
				sum_XY=np.sum(x_template*x_data)
				sum_X,sum_Y=np.sum(x_template),np.sum(x_data)
				cov=(m*sum_XY-sum_X*sum_Y)/(m*(m-1))

				X_square=np.sum(x_template**2)
				Y_square=np.sum(x_data**2)
				deno1=m*Y_square-np.sum(x_data)**2
				deno2=(m*X_square-np.sum(x_template)**2)
				vary=np.sqrt(deno2*deno1)/(m*(m-1))
				r_lag=cov/vary
				lag_arr.append(lag)
				r_arr.append(r_lag)





	return np.asarray(lag_arr),np.asarray(r_arr)	





