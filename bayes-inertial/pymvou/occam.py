from __future__ import division
import sys
import pymvou
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec


#gs = gridspec.GridSpec(nrows=1,ncols=2,height_ratios = [1,1],width_ratios=[1,1],wspace =0.38,hspace = 0.33)
#gs.update(top = 0.97, bottom = 0.12, left = 0.19, right = 0.90)
gs = gridspec.GridSpec(nrows=3,ncols=2,wspace =0.38,hspace = 0.33)

with sns.axes_style("white"):
    ax1 = plt.subplot(gs[0,0])
    ax1.axes.get_yaxis().set_ticks([])
with sns.axes_style("white"):
    ax2 = plt.subplot(gs[0,1])
    ax2.axes.get_yaxis().set_ticks([])
with sns.axes_style("white"):
    ax3 = plt.subplot(gs[1,0])
    ax3.axes.get_yaxis().set_ticks([])
with sns.axes_style("white"):
    ax4 = plt.subplot(gs[1,1])
    ax4.axes.get_yaxis().set_ticks([])
with sns.axes_style("white"):
    ax5 = plt.subplot(gs[2,0])
    ax5.axes.get_yaxis().set_ticks([])
with sns.axes_style("white"):
    ax6 = plt.subplot(gs[2,1])
    ax6.axes.get_yaxis().set_ticks([])

ax1 = plt.subplot(gs[0,0])
ind = (0.2,0.8)
width = 0.4

theta1 = 0.5
sigma1 = 0.1
x1, v1  = pymvou.ppou(theta1,sigma1)
lpo1 = pymvou.occamfactor(x1,v1)

ax1.bar(ind[0], lpo1[0], width, color = 'r',alpha = 0.7)
ax1.bar(ind[1], lpo1[1], width, color = 'c')
plt.xlim([0,1.4])
###
plt.ylim([0,200000])
###
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
plt.ylabel('logP(Evidence)')
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)
plt.title(r'$Overdamped$', fontsize = 30)

ax2 = plt.subplot(gs[0,1])
ind = (0.2,0.8)
width = 0.4

simparams2 = pymvou.pconverter(32768,1/65536,1e-12,225e-6,10)
x2, v2 = pymvou.simxv(simparams2, pymvou.sigmamatrix(simparams2), pymvou.explambda(simparams2))
lpo2 = pymvou.occamfactor(x2,v2)

ax2.bar(ind[0], lpo2[0], width, color = 'r', alpha = 0.7)
ax2.bar(ind[1], lpo2[1], width, color = 'c')
plt.xlim([0,1.4])
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
###
plt.ylim([0,2000000])
###
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)
plt.title(r'$Underdamped$', fontsize = 30)

ax3 = plt.subplot(gs[1,0])
ind = (0.2,0.8)
width = 0.4

theta3 = 0.7
sigma3 = 0.1
x3, v3  = pymvou.ppou(theta3,sigma3)
lpo3 = pymvou.occamfactor(x3,v3)

ax3.bar(ind[0], lpo3[0], width, color = 'r', alpha = 0.7)
ax3.bar(ind[1], lpo3[1], width, color = 'c')
plt.xlim([0,1.4])
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
###
plt.ylim([0,250000])
###
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)


ax4 = plt.subplot(gs[1,1])
ind = (0.2,0.8)
width = 0.4

simparams4 = pymvou.pconverter(32768,1/65536,1.5e-12,250e-6,20)
x4, v4 = pymvou.simxv(simparams4, pymvou.sigmamatrix(simparams4), pymvou.explambda(simparams4))
lpo4 = pymvou.occamfactor(x4,v4)

ax4.bar(ind[0], lpo4[0], width, color = 'r', alpha = 0.7)
ax4.bar(ind[1], lpo4[1], width, color = 'c')
plt.xlim([0,1.4])
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
###
plt.ylim([0,2000000])
###
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)


ax5 = plt.subplot(gs[2,0])
ind = (0.2,0.8)
width = 0.4

theta5 = 0.6
sigma5 = 0.1
x5, v5  = pymvou.ppou(theta5,sigma5)
lpo5 = pymvou.occamfactor(x5,v5)

ax5.bar(ind[0], lpo5[0], width, color = 'r', alpha = 0.7)
ax5.bar(ind[1], lpo5[1], width, color = 'c')
plt.xlim([0,1.4])
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
###
plt.ylim([0,250000])
###
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)


ax6 = plt.subplot(gs[2,1])
ind = (0.2,0.8)
width = 0.4

simparams6 = pymvou.pconverter(32768,1/65536,2e-12,300e-6,30)
x6, v6 = pymvou.simxv(simparams6, pymvou.sigmamatrix(simparams6), pymvou.explambda(simparams6))
lpo6 = pymvou.occamfactor(x6,v6)

ax6.bar(ind[0], lpo6[0], width, color = 'r', alpha = 0.7)
ax6.bar(ind[1], lpo6[1], width, color = 'c')
plt.xlim([0,1.4])
plt.xticks([0.4,1],[r'$H_{OD}$', r'$H_{UD}$'], fontsize = 30)
###
plt.ylim([0,2000000])
###
plt.ylabel(r'$lnP(D|H)$', fontsize = 30)


"""
def occamfactor(x,v):
	k_B = 1.38064881313131e-23  # Boltzmann constant (Newton metre/Kelvin)
	T = 275 #Kelvin
	N = len(x)
	k = (k_B*T)/np.var(x)
	m = (k_B*T)/np.var(v)
	log_ud = N*np.log(k/(2*np.pi*k_B*T))+N*np.log(m/(2*np.pi*k_B*T))-(2*N)+np.log((4*np.pi*m*k)/N)
	log_od = N*np.log(k/(2*np.pi*k_B*T))-(1*N)+np.log((2*np.pi*k)/np.sqrt(N))
	return log_ud, log_od

#x = np.fromfile('position.txt',dtype=float,count=-1,sep = " ")
#v = np.fromfile('velocity.txt',dtype = float, count=-1, sep = " ")
print occamfactor(np.fromfile('position.txt',dtype=float,count=-1,sep = " "),np.fromfile('velocity.txt',dtype = float, count=-1, sep = " "))
"""
plt.show()