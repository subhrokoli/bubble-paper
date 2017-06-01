from __future__ import division
import pymvou
import numpy as np
import matplotlib.pyplot as plt

#gamma_crit is the critical damping which is given by sqrt(4*m*k)
#for the simulations, gamma is taken as
#gamma = (gamma_crit/gamma_factor)
#this is how the value of gamma is tuned through the parameter gamma_factor

N = 32768 #Number of data points
f_sam = 65536 #Hz
dt = 1/f_sam
m = 1e-12 #1 nanogram
T = 275 #Kelvin
k = 225e-6 #kg/s2
gamma_factor = 10

simparams = pymvou.pconverter(N,dt,T,m,k,gamma_factor)

sigma_matrix = pymvou.sigmamatrix(simparams)
expM = pymvou.explambda(simparams)

x, v = pymvou.simxv(simparams, sigma_matrix, expM)
k_stat, m_stat, k_bayes, m_bayes, gamma_bayes = pymvou.estimator(x,v,T)
xo, vo = pymvou.ppou(0.5,0.1)

# k_stat is stationary estimate of trap stiffness
# m_stat is stationary estimate of mass
# k_bayes is Bayesian estimate of trap stiffness
# m_bayes is Bayesian estimate of mass
# gamma_bayes is Bayesian estimate of gamma


# x and v are position and velocity trajectories for underdamped data using parameters defined at the beginning of the code
# xo and vo are overdamped data generated using pyprocess
# the paramters of pymvou.ppou are (theta,sigma)
# where theta = lambda (dimensionless)
# and sigma = sqrt(2D), again dimensionless

#Insert code here to output any data into text files