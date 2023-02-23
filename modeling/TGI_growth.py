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
x = xg["x2"][2:-2]  # remove ghost cells
y = xg["x3"][2:-2]
z = xg["x1"][2:-2]

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
    
# attempt at computing approximate linear growth rate
dat=gemini3d.read.frame(direc,cfg["time"][0])   # get background state info directly from model initial conditions
ky=2*np.pi/1000                                  # visually the fastest growing mode is about 750 m
ne=dat["ne"][ialt,:,:]
n0=ne.mean()                                    # try grid-averaged values for these, use F-region ref. altitude
Ti=dat["Ti"][ialt,:,:]
ix3=int(np.floor(xg["lx"][2]/2))
Ti0=Ti.mean()
dx=x[ix2+1]-x[ix2-1]
dndx=(ne[ix2+1,ix3]-ne[ix2-1,ix3])/dx
dTdx=(Ti[ix2+1,ix3]-Ti[ix2-1,ix3])/dx
q=1.6e-19
B=45000e-9
kB=1.38e-23
gamma=ky/q/B*np.sqrt(kB*Ti0/n0*abs(dndx)*kB*abs(dTdx))     # corrected version of approximate growth from Keskinen, 2004

# define a reference exponential growth profile
tref=500                                         # reference time after initial noise has settled and growth has begun
itref=np.argmin(abs(tamp-tref))
nref=neamp[itref];
nlinear=nref*np.exp(np.array(gamma)*(np.array(tamp)-tref))

# plots
plt.subplots(2,1,dpi=250)
plt.subplot(2,1,1)
plt.pcolormesh(tamp[0:it],y,neslices[0:it,:].transpose())
plt.colorbar()
plt.xlabel("time (s)")
plt.ylabel("$\Delta n_e$ (m$^{-3}$)")
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
