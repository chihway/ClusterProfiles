
import numpy as np
import os 
import sys

# this is to combine all the JK files and output the mean and covariance

name = sys.argv[1]
Njk = int(sys.argv[2])
result_dir = sys.argv[3]

R = np.load(result_dir+name+'/Sigmag_0.npz')['R']
Sg = []

for i in range(Njk):
    mean = np.load(result_dir+name+'/Sigmag_'+str(i)+'.npz')['mean']
    Sg.append(mean)

Sg = np.array(Sg)
C = np.cov(Sg.T)*(len(Sg)-1)
Sg_mean = np.mean(Sg, axis=0)
Sg_sig = np.sum((Sg-Sg_mean)**2, axis=0)**0.5/len(Sg)**0.5*(len(Sg)-1)**0.5

np.savez('splashback_cov_'+str(name)+'.npz', cov=C,
         r_data=R, sg_mean=Sg_mean, sg_sig=Sg_sig)

