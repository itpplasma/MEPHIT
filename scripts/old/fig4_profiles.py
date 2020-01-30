# -*- coding: utf-8 -*-
"""
Created on Sat Aug 06 08:48:24 2016

@author: chris
"""
import numpy as np
import matplotlib.pyplot as plt

from exportfig import exportfig

plt.close('all')
plt.rc('font',**{'family':'serif','serif':['Times'],'size':'8'})
plt.rc('text', usetex=True)
plt.rc('figure', figsize=[6.4, 2.13])
plt.rc('lines', markersize=3)
colors = plt.rcParams['axes.color_cycle']

datafile1 = 'hpsi_vac.spectrum'
datafile2 = 'hpsi_relaxed_71.spectrum'

data1 = np.loadtxt(datafile1)
data2 = np.loadtxt(datafile2)
#data3 = np.loadtxt(datafile3)
r = np.sqrt((np.max(data1[:,0])-data1[:,0])/np.max(data1[:,0]))**2
#s3 = data3[:,0]
q = data1[:,1]

Nm1 = len(data1[0,2:])/2
bpsi1 = data1[:,2:(2+Nm1)] + 1j*data1[:,(2+Nm1):]
Nm2 = len(data2[0,2:])/2
bpsi2 = data2[:,2:(2+Nm2)] + 1j*data2[:,(2+Nm2):]

plt.figure()
for k in range(8):
    plt.subplot(2,4,k+1)
    #plt.title('m={}'.format(k))
    plt.plot(r, np.abs(bpsi1[:,(Nm1-1)/2+k]+bpsi2[:,(Nm1-1)/2+k]))
    plt.plot(r, np.abs(bpsi1[:,(Nm1-1)/2+k]),'r--')
    yrang = [0,250]
    if(k<4):
        yrang = [0,100]
        
    plt.plot(np.array([0.608,0.608])**2, yrang,'b',alpha=.5)#1.5
    plt.plot(np.array([0.760,0.760])**2, yrang,'b',alpha=.5)#2
    plt.plot(np.array([0.823,0.823])**2, yrang,'b',alpha=.5)#2.5
    plt.plot(np.array([0.861,0.861])**2, yrang,'b',alpha=.5)#3
    plt.plot(np.array([0.891,0.891])**2, yrang,'b',alpha=.5)#3.5
    plt.plot(np.array([0.918,0.918])**2, yrang,'b',alpha=.5)#4
    #plt.plot([0.250,0.250], yrang,'g-')
    plt.gca().text(0.1,.8*yrang[1],'m={}'.format(k))
    plt.gca().text(0.43,-0.33*yrang[1],r'$\rho_{\rm pol}^2$')
    plt.gca().text(-0.31,0.54*yrang[1],r'$\left\vert \delta B^{\psi}_{m}\right\vert$',
    rotation='vertical')
    #plt.legend(['m={}'.format(k)],)
    #if(k>0):
    #    plt.gca().xaxis.set_major_formatter(plt.NullFormatter())
    #    plt.gca().yaxis.set_major_formatter(plt.NullFormatter())
    #else:
    #plt.xlabel(r'$\rho_{\rm pol}$')
    #plt.ylabel(r'$B_r$')
    plt.ylim(yrang)
    #plt.grid(True)
    plt.tight_layout()
    
exportfig('fig4_profiles')