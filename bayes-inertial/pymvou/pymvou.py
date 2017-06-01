from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.optimize import fsolve
sys.path.append('/home/dipanjan/PyProcess/pyprocess')
import pyprocess as pp

#simulated data
k_B = 1.38064881313131e-23  # Boltzmann constant (Newton metre/Kelvin)
#T = 275 #Kelvin
"""
N = 32768 #Number of data points
m = 1e-12 #1 nanogram
T = 275 #Kelvin
k = 225e-6 #kg/s2
gamma_factor = 10
f_sam = 65536 #Hz
dt = 1/f_sam
gamma2_crit = 4*m*k
gamma_crit = np.sqrt(4*m*k)
gamma2 = gamma2_crit/(gamma_factor**2) #
gamma = gamma_crit/gamma_factor #
omega0 = np.sqrt(k/m)
tau = m/gamma
omega = np.sqrt((omega0**2)-(1/(4*tau*tau)))
k_B = 1.38064881313131e-23  # Boltzmann constant (Newton metre/Kelvin)
D = k_B*T/gamma
simparams = (N, dt, omega, omega0, tau, D)
"""
def pconverter(N, dt, T, mass, spring, gamma_f):
    m = mass
    k = spring
    gamma_factor = gamma_f

    gamma2_crit = 4*m*k
    gamma_crit = np.sqrt(4*m*k)
    gamma2 = gamma2_crit/(gamma_factor**2) #
    gamma = gamma_crit/gamma_factor #
    omega0 = np.sqrt(k/m)
    tau = m/gamma
    omega = np.sqrt((omega0**2)-(1/(4*tau*tau)))
    D = k_B*T/gamma
    simparams = (N, dt, omega, omega0, tau, D)
    return N, dt, T, omega, omega0, tau, D



#Function to costruct sigma matrix
def sigmamatrix(simparams):
    N, dt, T, omega, omega0, tau, D = simparams
    ss1 = np.cos(2*omega*dt)-(2*omega*tau*np.sin(2*omega*dt))-(4*(omega0**2)*(tau**2))
    ss2 = np.cos(2*omega*dt)+(2*omega*tau*np.sin(2*omega*dt))-(4*(omega0**2)*(tau**2))

    sigma2_xx = (D/(4*(omega**2)*(omega0**2)*(tau**3)))*((4*(omega**2)*(tau**2))+(np.exp(-dt/tau))*(ss1))
    sigma2_vv = (D/(4*(omega**2)*(tau**3)))*((4*(omega**2)*(tau**2))+(np.exp(-dt/tau))*(ss2))
    sigma2_xv = (D/((omega**2)*(tau**2)))*(np.exp(-dt/tau)*np.sin(omega*dt)*np.sin(omega*dt))
    
    return sigma2_xx, sigma2_vv, sigma2_xv    

#sigma_matrix = sigmamatrix(simparams)

#Function to construct exponential matrix
def explambda(simparams):
    N, dt, T, omega, omega0, tau, D = simparams
    I = np.eye(2)
    J11 =(1/(2*omega*tau))
    J12 = (1/omega)
    J21 = -(omega0**2)/omega
    J22 = -J11
    J = np.matrix([[J11,J12],[J21,J22]])

    return np.exp(-dt/(2*tau))*((np.cos(omega*dt)*I)+(np.sin(omega*dt)*J))
    
#expM = explambda(simparams)

def simxv(simparams, sigma_matrix, expM):
    N, dt, T, omega, omega0, tau, D = simparams
    x = np.zeros([N,1])
    v = np.zeros([N,1])
    
    sigma2_xx, sigma2_vv, sigma2_xv = sigma_matrix

    for j in np.arange(0,N-1):
        oldvec = np.array([x[j],v[j]])
        randgauss = np.random.randn(2,1)
        delx = np.sqrt(sigma2_xx)*randgauss[0]
        delv = (sigma2_xv/(np.sqrt(sigma2_xx)))*randgauss[0]+(np.sqrt(sigma2_vv - ((sigma2_xv**2)/(sigma2_xx))))*randgauss[1]
        delvec = np.array([delx,delv])
        updatevec = np.dot(expM,oldvec)+delvec
        x[j+1] = updatevec[0]
        v[j+1] = updatevec[1]
    return x,v

#Estimates
def estimator(x,v,T):
    k_B = 1.38064881313131e-23  # Boltzmann constant (Newton metre/Kelvin)
    N = np.prod(x.shape) #length of dataset   
    
    #Stationary estimates
    k_stat = k_B*T/(np.var(x))
    m_stat = k_B*T/(np.var(v))
    
    #Bayesian estimates
    #matrix sufficient statistics
    T1_11 = np.dot(np.transpose(x[1:N]),x[1:N])[0,0]
    T1_12 = np.dot(np.transpose(x[1:N]),v[1:N])[0,0]
    T1_21 = T1_12
    T1_22 = np.dot(np.transpose(v[1:N]),v[1:N])[0,0]

    T2_11 = np.dot(np.transpose(x[0:N-1]),x[1:N])[0,0]
    T2_12 = np.dot(np.transpose(x[1:N]),v[0:N-1])[0,0]
    T2_21 = np.dot(np.transpose(x[0:N-1]),v[1:N])[0,0]
    T2_22 = np.dot(np.transpose(v[0:N-1]),v[1:N])[0,0]

    T3_11 = np.dot(np.transpose(x[0:N-1]),x[0:N-1])[0,0]
    T3_12 = np.dot(np.transpose(x[0:N-1]),v[0:N-1])[0,0]
    T3_21 = T3_12
    T3_22 = np.dot(np.transpose(v[0:N-1]),v[0:N-1])[0,0]


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

def occamfactor(x,v,T):
    N = len(x)
    k = (k_B*T)/np.var(x)
    m = (k_B*T)/np.var(v)
    log_ud = N*np.log(k/(2*np.pi*k_B*T))+N*np.log(m/(2*np.pi*k_B*T))-(2*N)+np.log((4*np.pi*m*k)/N)
    log_od = N*np.log(k/(2*np.pi*k_B*T))-(1*N)+np.log((2*np.pi*k)/np.sqrt(N))
    return log_od, log_ud

def ppou(theta,sigma):
    #theta mean reversion rate
    mu = 0.0     # secular trend
    #sigma volatility
    r = 0.1             # dimensionless sampling rate
    M = 10000           # dimensionless sampling duration
    #M = 0.1*2**16
    dt = r * theta      # sampling rate
    tf = M * theta      # sampling duration
    N = int(M/r)        # number of sample points
    t = dt*np.arange(N) # sampling times
    oup = pp.OU_process(theta=theta, mu=mu, sigma=sigma, startPosition=0)
    x = oup.sample_path(t)
    x = x.flatten()
    oup_data = { 'N':N, 'dt':dt, 'x':x}
    v = np.zeros([N,1])
    for index in np.arange(1,N):
        v[index] = (1/dt)*(x[index]-x[index-1])
    return x, v

#x, v = simxv(simparams, sigma_matrix, expM)
#x, v = simxv(simparams, sigmamatrix(simparams), expM)
"""
x.tofile('position.txt', sep=" ", format="%s")
v.tofile('velocity.txt', sep=" ", format="%s")
"""