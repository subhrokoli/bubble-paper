from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.optimize import fsolve

#simulated data
N = 32768 #Number of data points
m = 2e-12 #1 nanogram
T = 300 #Kelvin
k = 300e-6 #kg/s2
gamma_factor = 30
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
simparams = (N, dt, omega, omega0, tau)

#Function to costruct sigma matrix
def sigmamatrix(simparams):
    simparams = (N, dt, omega, omega0, tau)
    ss1 = np.cos(2*omega*dt)-(2*omega*tau*np.sin(2*omega*dt))-(4*(omega0**2)*(tau**2))
    ss2 = np.cos(2*omega*dt)+(2*omega*tau*np.sin(2*omega*dt))-(4*(omega0**2)*(tau**2))

    sigma2_xx = (D/(4*(omega**2)*(omega0**2)*(tau**3)))*((4*(omega**2)*(tau**2))+(np.exp(-dt/tau))*(ss1))
    sigma2_vv = (D/(4*(omega**2)*(tau**3)))*((4*(omega**2)*(tau**2))+(np.exp(-dt/tau))*(ss2))
    sigma2_xv = (D/((omega**2)*(tau**2)))*(np.exp(-dt/tau)*np.sin(omega*dt)*np.sin(omega*dt))
    
    return sigma2_xx, sigma2_vv, sigma2_xv    

sigma_matrix = sigmamatrix(simparams)

#Function to construct exponential matrix
def explambda(simparams):
    N, dt, omega, omega0, tau = simparams
    I = np.eye(2)
    J11 =(1/(2*omega*tau))
    J12 = (1/omega)
    J21 = -(omega0**2)/omega
    J22 = -J11
    J = np.matrix([[J11,J12],[J21,J22]])

    return np.exp(-dt/(2*tau))*((np.cos(omega*dt)*I)+(np.sin(omega*dt)*J))
    
expM = explambda(simparams)

def simxv(simparams, sigma_matrix, expM):
    N, dt, omega, omega0, tau = simparams
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

x, v = simxv(simparams, sigma_matrix, expM)

x.tofile('position.txt', sep=" ", format="%s")
v.tofile('velocity.txt', sep=" ", format="%s")