#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:37:17 2023

This script will read in a GEMINI TGI instability simulation and attempt to plot the
growth of density irregularities and compare it against linear theory.  

dependencies:  pygemini

@author: zettergm
"""

# imports
import gemini3d.read
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
#from plotGDI_tools import padstr

# set some font sizes
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# location of simulation output
home=os.path.expanduser("~")
direc = home+"/simulations/TGI_MR_veryhighne_norecomb/"
# plotdir=direc+"/customplots/"
# if not os.path.isdir(plotdir):
#     os.mkdir(plotdir)
parmlbl="ne"

# config and grid data
print("...Loading config and grid...")
cfg = gemini3d.read.config(direc)
xg = gemini3d.read.grid(direc)
x = np.array(xg["x2"][2:-2])  # remove ghost cells
y = np.array(xg["x3"][2:-2])
z = np.array(xg["x1"][2:-2])

# reference altitude
altref = 300e3
ialt = np.argmin(abs(z - altref), axis=0)

# input electric field and drifts
t0 = cfg["time"][0]

# load data from a specified set of time indices
its=range(0,len(cfg["time"]),1)
ix2=int(np.floor(xg["lx"][1]/2))
# plt.ioff()
# plt.figure(dpi=150)
neamp=np.zeros(len(its))
tamp=np.zeros(len(its))
neslices=np.zeros((len(its),xg["lx"][2]))
for it in its:
    print("Loading:  ",cfg["time"][it])
    dat=gemini3d.read.frame(direc,cfg["time"][it])
    ne=np.array(dat[parmlbl])
    neslices[it,:]=ne[ialt,ix2,:]-ne[ialt,ix2,:].mean()    # subtract out mean so we have only fluctuations
    tamp[it]=(cfg["time"][it]-t0).total_seconds()
    neamp[it]=np.std(neslices[it,:]);
N=it     # save successful number of time samples loaded
    
# attempt at computing approximate linear growth rate
dat=gemini3d.read.frame(direc,cfg["time"][0])   # get background state info directly from model initial conditions
ky=2*np.pi/1000                                  # visually the fastest growing mode is about 750 m
ne=np.array(dat["ne"][ialt,:,:])
n0=ne.mean()                                    # try grid-averaged values for these, use F-region ref. altitude
Ti=np.array(dat["Ti"][ialt,:,:])
ix3=int(np.floor(xg["lx"][2]/2))
Ti0=Ti.mean()
dx=x[ix2+1]-x[ix2-1]
dndx=(ne[ix2+1,ix3]-ne[ix2-1,ix3])/dx
dTdx=(Ti[ix2+1,ix3]-Ti[ix2-1,ix3])/dx
q=1.6e-19
B=45000e-9
kB=1.38e-23
gamma=ky/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))     # corrected version of approximate growth from Keskinen, 2004

# compute some derived parameters like effective scale lengths for density and temperature variation
LT=1/(1/Ti0*dTdx)
Ln=1/(1/n0*dndx)
print("Effective scale lengths at sim. beginning:  ",LT,Ln)

# define a reference exponential growth profile
tref=500                                         # reference time after initial noise has settled and growth has begun
itref=np.argmin(abs(tamp-tref))
nref=neamp[itref];
nlinear=nref*np.exp(gamma*(tamp-tref))

# plots of spatial domain fluctuations
plt.subplots(2,1,dpi=150)
plt.subplot(2,1,1)
plt.pcolormesh(tamp[0:it],y,neslices[0:it,:].transpose()/n0)
#plt.pcolormesh(tamp[0:it],y,neslices[0:it,:].transpose())
plt.colorbar()
plt.xlabel("time (s)")
plt.title("$\Delta n_e / n_e$")
plt.ylabel("dist. ortho. to gradients (m)")
#plt.ylabel("$\Delta n_e (m$^{-3}$)")
plt.subplot(2,1,2)
plt.plot(tamp[0:it],np.log10(neamp[0:it]))
plt.xlabel("time (s)")
plt.ylabel("$log_{10}~~\sigma_{\Delta n_e}$ (m$^{-3}$)")
axes=plt.gca()
ylims=axes.get_ylim()
plt.plot(tamp[0:it],np.log10(nlinear[0:it]))
plt.ylim(ylims)     # reset axes to something not crazy
plt.legend(("model output","linear theory from K04"))
plt.show()

# Spectral analysis of fluctuations
dy=y[1]-y[0]        # spacing
ly=int(len(y))
xf=scipy.fft.fftfreq(ly,dy)   # frequency axis for ffts
nf=np.empty((N,ly),dtype=complex)
nfpwr=np.zeros((N,ly))
for it in range(N):
    nf[it,:]=scipy.fft.fft(neslices[it,:])
    nfpwr[it,:]=nf[it,:]*np.conjugate(nf[it,:])

# only use positive frequencies (since real-valued signal negatives are linearly dependent)
xfplus=xf[0:ly//2]
nfpwrplus=nfpwr[0:N,0:ly//2]
totalpwrplus=np.sum(nfpwrplus,1)

# Select a few reference frequencies to track
ifreq1=np.argmin(abs(xfplus-0.00135))    # fastest growing mode
ifreq2=np.argmin(abs(xfplus-0.0006))     # longer wavelength feature? sub-harmonic???
ifreq3=np.argmin(abs(xfplus-0.0002))     # larger feature appearing near end?
ifreq4=np.argmin(abs(xfplus-0.0018))     # features that loses energy in the middle

# plots of frequency domain fluctuations
plt.subplots(2,1,dpi=150)
plt.subplot(2,1,1)
#plt.pcolormesh(tamp[0:N],xf[0:ly//2],np.abs(nf[0:N,0:ly//2].transpose()))
plt.pcolormesh(tamp[0:N],xfplus,np.log10(nfpwrplus.transpose()))
plt.xlabel('time (s)')
plt.ylabel('wavenumber (m$^{-1}$)')
plt.title("$\Delta n_e$ power spectrum (log$_{10})$")
plt.colorbar()
plt.ylim(0,0.005)
plt.clim(21,23)
plt.subplot(2,1,2)
plt.plot(tamp[0:N],np.log10(nfpwrplus[:,ifreq1]))
plt.plot(tamp[0:N],np.log10(nfpwrplus[:,ifreq2]))
plt.plot(tamp[0:N],np.log10(nfpwrplus[:,ifreq3]))
plt.plot(tamp[0:N],np.log10(nfpwrplus[:,ifreq4]))
plt.plot(tamp[0:N],np.log10(totalpwrplus))
plt.legend((str(1/xf[ifreq1]), str(1/xf[ifreq2]), str(1/xf[ifreq3]), str(1/xf[ifreq4]), "total" ))
plt.xlabel("time (s)")
plt.ylabel("$\Delta n_e$ power spectrum (log$_{10})$")
plt.show()

# Find reference density fluctuation amplitudes at these times, (note sqrt to go back to density units)
n1=abs(nf[itref,ifreq1])
n2=abs(nf[itref,ifreq2])
n3=abs(nf[itref,ifreq3])
n4=abs(nf[itref,ifreq4])

# Evaluate the growth rate at various wavenumbers
k1=2*np.pi*xf[ifreq1]
k2=2*np.pi*xf[ifreq2]
k3=2*np.pi*xf[ifreq3]
k4=2*np.pi*xf[ifreq4]

# Now distinct growth rates
gamma1=k1/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))
gamma2=k2/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))
gamma3=k3/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))
gamma4=k4/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))

# Now exponential profiles 
nlinear1=n1*np.exp(gamma1*(tamp-tref))
nlinear2=n2*np.exp(gamma2*(tamp-tref))
nlinear3=n3*np.exp(gamma3*(tamp-tref))
nlinear4=n4*np.exp(gamma4*(tamp-tref))

# Now plot against unstable growth modes
plt.subplots(4,1, dpi=150)

plt.subplot(4,1,1)
plt.plot(tamp[0:N],np.log10(abs(nf[0:N,ifreq1])))
axes=plt.gca()
ylims=axes.get_ylim()
plt.plot(tamp,np.log10(nlinear1))
plt.ylim(ylims)
plt.xlabel("time (s)")
plt.ylabel("Fourier coeff.")
plt.title(str(2*np.pi/k1))

plt.subplot(4,1,2)
plt.plot(tamp[0:N],np.log10(abs(nf[0:N,ifreq2])))
axes=plt.gca()
ylims=axes.get_ylim()
plt.plot(tamp,np.log10(nlinear2))
plt.ylim(ylims)
plt.xlabel("time (s)")
plt.ylabel("Fourier coeff.")
plt.title(str(2*np.pi/k2))

plt.subplot(4,1,3)
plt.plot(tamp[0:N],np.log10(abs(nf[0:N,ifreq3])))
axes=plt.gca()
ylims=axes.get_ylim()
plt.plot(tamp,np.log10(nlinear3))
plt.ylim(ylims)
plt.xlabel("time (s)")
plt.ylabel("Fourier coeff.")
plt.title(str(2*np.pi/k3))

plt.subplot(4,1,4)
plt.plot(tamp[0:N],np.log10(abs(nf[0:N,ifreq4])))
axes=plt.gca()
ylims=axes.get_ylim()
plt.plot(tamp,np.log10(nlinear4))
plt.ylim(ylims)
plt.xlabel("time (s)")
plt.ylabel("Fourier coeff.")
plt.title(str(2*np.pi/k4))
plt.show()


