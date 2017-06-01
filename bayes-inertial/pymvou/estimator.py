from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.optimize import fsolve

T = 300 #Kelvin
x = np.fromfile('position.txt',dtype=float,count=-1,sep = " ")
v = np.fromfile('velocity.txt',dtype = float, count=-1, sep = " ")

#Estimates
def estimator(x,v,T):
    k_B = 1.38064881313131e-23  # Boltzmann constant (Newton metre/Kelvin)
    N = np.prod(x.shape) #length of dataset   
    
    #Stationary estimates
    k_stat = k_B*T/(np.var(x))
    m_stat = k_B*T/(np.var(v))
    
    #Bayesian estimates
    #matrix sufficient statistics
    T1_11 = np.dot(np.transpose(x[1:N]),x[1:N])
    T1_12 = np.dot(np.transpose(x[1:N]),v[1:N])
    T1_21 = T1_12
    T1_22 = np.dot(np.transpose(v[1:N]),v[1:N])

    T2_11 = np.dot(np.transpose(x[0:N-1]),x[1:N])
    T2_12 = np.dot(np.transpose(x[1:N]),v[0:N-1])
    T2_21 = np.dot(np.transpose(x[0:N-1]),v[1:N])
    T2_22 = np.dot(np.transpose(v[0:N-1]),v[1:N])

    T3_11 = np.dot(np.transpose(x[0:N-1]),x[0:N-1])
    T3_12 = np.dot(np.transpose(x[0:N-1]),v[0:N-1])
    T3_21 = T3_12
    T3_22 = np.dot(np.transpose(v[0:N-1]),v[0:N-1])


    T1 = np.matrix([[T1_11,T1_12],[T1_21,T1_22]])
    T2 = np.matrix([[T2_11,T2_12],[T2_21,T2_22]])
    T3 = np.matrix([[T3_11,T3_12],[T3_21,T3_22]])

    invT3 = np.linalg.inv(T3)
    
    Sigma_est = (1/N)*(T1 - (T2*invT3*np.transpose(T2)))
    expMest = T2*invT3

    coeffmat = np.matrix([[1-(expMest[0,0]**2),-expMest[0,1]**2],[-expMest[1,0]**2,1-(expMest[1,1]**2)]])
    invcoeffmat = np.linalg.inv(coeffmat)
    sigvec = np.array([Sigma_est[0,0],Sigma_est[1,1]])
    cvec = np.dot(invcoeffmat,sigvec)
    k_bayes = (k_B*T)/cvec[0,0]
    m_bayes = (k_B*T)/cvec[0,1]
    
    gamma_bayes = (Sigma_est[0,1]*(m_bayes**2))/(k_B*T*(expMest[0,1]**2))
    
    return k_stat, m_stat, k_bayes, m_bayes, gamma_bayes

k_stat, m_stat, k_bayes, m_bayes, gamma_bayes = estimator(x,v,T)
print k_stat
print m_stat
print k_bayes
print m_bayes
print gamma_bayes