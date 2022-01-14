# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:26:39 2021

@author: sid
"""
from matplotlib import pyplot as plt
import pandas
import os
import numpy as np
import plotly.express as px
import plotly.offline as po
import plotly.graph_objs as go
#Use of fabio for IO and plotting
from silx.resources import ExternalResources
from scipy import ndimage, signal
import math
plt.ion()
#Detector calibration and setup
import pyFAI, pyFAI.detectors, fabio
print("Using pyFAI version", pyFAI.version)
import pyFAI.distortion as dis
from pyFAI.gui import jupyter
from pyFAI.calibrant import get_calibrant
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import time
import lmfit
import pandas as pd

#openfile for data #Only RAW DATA has radial averaging, as urwarped is already averaged
fig_z = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\Aerogel5_2_60s_107cm_01_unwarped.gfrm") #data has additional commands like shape
img_z = fig_z.data
fig_x = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\Aerogel1_1_60s_107cm_01_unwarped.gfrm")
img_x = fig_x.data
Mask_dat = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\air_60s_107cm_01_unwarped.gfrm")
msk = Mask_dat.data #Mask_correction
#Flats and background
flat_xz = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\2D-A5.2_01_001.gfrm") #Masks
flat = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\2D-A5.2_01_000.gfrm") #Marks
# fig_plt, ax = plt.subplots()
#1 = 0%
#2 = 85%
#3 = 40%
#4 = 70%
#5 = 98%
#6 = 95%
#for headers -> fig_z.header()
#Detector and measurement parameters
wl = 1.5418e-10 #nm dimension of X-rays
cal = get_calibrant("AgBh") #Silver behanate sample
cal.wavelength=wl
start_time = time.time()
print("PyFAI version", pyFAI.version)
Vantec = pyFAI.detectors.Detector(68e-6, 68e-6)#pixel size
Vantec.max_shape=(2048,2048)#image shape
ai = AzimuthalIntegrator(dist=1.07050, detector=Vantec, wavelength=wl)#initialization of arbitrary detector with given dimensions

#Masks and darks
Ai_mask = ai.create_mask(msk)
ai.mask = Ai_mask
#image center calculation
cent = msk.T
x_cent = np.zeros(len(cent))
x_holder = np.zeros(len(cent))
y_cent = np.zeros(len(cent))
y_holder = np.zeros(len(cent))
#image center calculation
for i in range(len(cent)):
    for j in range(len(cent)):
        x_holder[j] = cent[i][j] #running X center intensity loop
        y_holder[j] = cent[j][i] #running Y center intensity loop
    x_cent[i] = x_holder.sum()
    y_cent[i] = y_holder.sum()
x_c=y_c = 0
for i in range(len(cent)):
    ctr_x = (x_cent[i]*i)
    ctr_y = (y_cent[i]*i)
    x_c+=ctr_x
    y_c+=ctr_y

xx_ctr = x_c/x_cent.sum()
yy_ctr = y_c/y_cent.sum() #weighted average for center position

#finding beamcenter with PONI=Point of normal incedence
p1 = 68e-6 * 2048/2
ai.poni1 = p1 - 0.00017
p2 = 68e-6 * 2048/2
ai.poni2 = p2
print(ai)
fake = cal.fake_calibration_image(ai)
#detector setup complete
#Fixing the peak processing with a centering kernel
size = 11 #Odd of course
center = (size-1)//2
y, x = np.ogrid[-center:center+1,-center:center+1]
r2 = x*x + y*y
kernel = (r2<=(center+0.5)**2).astype(float)
kernel /= kernel.sum()
fig_fix,ax_fix = plt.subplots()
# ax_fix.imshow(kernel, interpolation="nearest", origin="lower")
cnv = signal.convolve2d(img_z, kernel, mode="same")
#convolution complete_but errors are still present

#CORRECTION FACTORS:
npt = 1000
kwarg = {"npt":npt,
          "correctSolidAngle":True,
          "polarization_factor":None,
          "safe":False}
omega  = ai.solidAngleArray(Vantec.shape, absolute=True)
flat = np.ones(Vantec.shape)
res_flat1 = ai.integrate1d(flat, 1000)
res_flat2 = ai.integrate2d(flat, 1000)
# crv = jupyter.plot1d(res_flat1)
# crv.axes.set_xlim(-1,15)
# crv.axes.set_ylim(0.9,1.1)
# crv2 = jupyter.plot2d(res_flat2)
#distortion correction
distort = dis.Distortion(detector=Vantec, shape=Vantec.shape, resize = False, empty=0,mask=msk,method = 'lut')
cor_img = distort.correct_ng(img_z, solidangle = ai.solidAngleArray)

#CORRECTIONS end:
plt.rcParams['figure.dpi'] = 600 #inline plot dpi setting
plt.rcParams["figure.figsize"] = (10,10)
Desc = 1500 #Descretizer for 2D unwrapping of the scattering
ai.setSPD(SampleDistance = 1.070500031,  Center_1=1024.6-2, Center_2 = 1026.5-0.5) #maybe removed, A2_2
#XZ plane plot
res_z = ai.integrate1d(img_z, 2000, unit="q_nm^-1",filename= "integrated.dat", radial_range = [0.1,3.0], correctSolidAngle=False, mask = msk)
rad_z = ai.integrate_radial(img_z, 2000, radial_range= [0.2,2.5],unit = "chi_deg", correctSolidAngle=True, mask = msk, npt_rad = 200)
jupyter.plot1d(res_z, label = "Compression axis")
jupyter.plot1d(rad_z)

Imat, Qv, tth = ai.integrate2d(img_z, 500, npt_azim=180, radial_range=[0.1,3.0], correctSolidAngle=False,mask=msk)

# jupyter.plot2d(d_x)
new= ndimage.rotate(img_z, 47, reshape=False)
new_1 = new[900:1150, 1024:1150]
plo = plt.imshow(new_1, cmap = 'magma', vmin=0, vmax=1000)

#2D Scatter plot without imshow
# jupyter.display(proc_img_z)#2d scatter display
# jupyter.display(proc_img_x)

plt.figure(figsize=(16,9))
plt.axis("off")
z = plt.imshow(new_1, cmap = 'gnuplot2', vmin=0, vmax = 1000)
plt.colorbar() # Show color bar of above image
plt.show()

def lmfit_Int_half(test, test1):
    gauss1 = lmfit.models.LorentzianModel(prefix='g1_')
    pars=gauss1.guess(test, x=test1)
    pars['g1_center'].set(value=5, min=-180, max=180)
    pars['g1_sigma'].set(value=200, min=-100)
    pars['g1_amplitude'].set(value=100, min=1)
    init = gauss1.eval(pars, x=test1)
    out = gauss1.fit(test, pars, x=test1)
    fwhm = 2*out.best_values['g1_sigma']*np.sqrt(2*np.log(2))
    # gauss2 = lmfit.models.LorentzianModel(prefix='g2_')
    # pars2=gauss2.guess(firstpeak, x=firsttth)
    # pars2['g2_center'].set(value=45, min=35, max=65)
    # pars2['g2_sigma'].set(value=120, min=40)
    # pars2['g2_amplitude'].set(value=100, min=10)
    # out_2 = gauss2.fit(firstpeak, pars2, firsttth)
    pie = (180-fwhm)/180
    # print(out.fit_report(min_correl=0.01))
    return out.best_fit, fwhm

def lmfit_Int(test, test1, peak1, tth1):
    gauss1 = lmfit.models.LorentzianModel(prefix='g1_')
    pars=gauss1.guess(test, x=test1)
    pars['g1_center'].set(value=5, min=-180, max=180)
    pars['g1_sigma'].set(value=500, min=-100)
    pars['g1_amplitude'].set(value=1000, min=1)
    init = gauss1.eval(pars, x=test1)
    out = gauss1.fit(test, pars, x=test1)
    fwhm = 2*out.best_values['g1_sigma']*np.sqrt(2*np.log(2))
    gauss2 = lmfit.models.LorentzianModel(prefix='g2_')
    pars2=gauss2.guess(peak1, x=tth1)
    pars2['g2_center'].set(value=5, min=-180, max=180)
    pars2['g2_sigma'].set(value=200, min=-100)
    pars2['g2_amplitude'].set(value=100, min=1)
    init_2 = gauss2.eval(pars2, x=tth1)
    out_2 = gauss2.fit(peak1, pars2, x=tth1)
    fwhm2 = 2*out_2.best_values['g2_sigma']*np.sqrt(2*np.log(2))
    pie = (180-fwhm)/180
    # print(out.fit_report(min_correl=0.01))
    return out.best_fit, out_2.best_fit, fwhm,fwhm2
# plt.plot(tth[0:65],a[0])
# plt.plot(tth[65:170], b[0])

ctr=ct=c=0
ctri=np.zeros(18)
holder= np.zeros(18)
holder_1 = np.zeros(18)
hold_q=np.zeros(18)
start = 0.3 #0.25--0.1, step 0.07
for i in range(len(tth)):
    if (Qv[i]>start and Qv[i]<1.1):
        ctri[ct]=i
        hold_q[ct] = Qv[i]
        start+=0.06
        ct+=1
       # break
Icum=np.zeros((len(tth),ct))
#     plt.plot(tth[170:290], Imat[170:290:,int(ctri[j])])
for i in range(ct):
    for j in range(len(tth)):
        Icum[j][c]+=Imat[j][int(ctri[c])]
    c+=1    
# ct=ct-1

n=np.zeros((len(tth), ct))
nn =np.zeros((len(tth), ct))

for i in range(len(tth)):
    hol=0
    for j in range(ct):
        hol += Icum[i][j]
        n[i][j]=hol  #actual cumulative

for i in range(ct):
    nn[0:170,i] = n[0:170,i]-n[0:170,i].min()
    
# for i in range(ct):
#     plt.plot(tth[90:145], Icum[90:145,i]) #120° profile for lorentzian fit
#     plt.plot(tth[90:145], lmfit_Int(Icum[90:145,i], tth[90:145])[0])
#     holder[i] = lmfit_Int(Icum[90:145,i], tth[90:145])[1]
    
# for i in range(ct):
#     plt.scatter(tth[0:170], n[0:170,i]) #120° profile for lorentzian fit
#     # plt.plot(tth[0:170], np.hstack((lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[0], lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[1])))
#     plt.plot(tth[0:60], np.hstack((lmfit_Int(n[0:60,i], tth[0:60], n[60:170,i], tth[60:170])[0])))
#     plt.plot(tth[75:155], np.hstack((lmfit_Int(n[0:65,i], tth[0:65], n[75:155,i], tth[75:155])[1])))
#     holder[i] = lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[2]
#     holder_1[i] = lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[3]
fig, ax = plt.subplots()

for i in range(ct):
    plt.scatter(tth[0:170], nn[0:170,i]) #120° profile for lorentzian fit
    # plt.plot(tth[0:170], np.hstack((lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[0], lmfit_Int(n[0:65,i], tth[0:65], n[65:170,i], tth[65:170])[1])))
    plt.plot(tth[0:60], np.hstack((lmfit_Int(nn[0:60,i], tth[0:60], nn[60:170,i], tth[60:170])[0])))
    plt.plot(tth[75:155], np.hstack((lmfit_Int(nn[0:65,i], tth[0:65], nn[75:155,i], tth[75:155])[1])))
    holder[i] = lmfit_Int(nn[0:65,i], tth[0:65], nn[65:170,i], tth[65:170])[2]
    holder_1[i] = lmfit_Int(nn[0:65,i], tth[0:65], nn[65:170,i], tth[65:170])[3]
    plt.xlim([-180, 180])
    plt.ylim([0, 500])
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tick_params(axis='both', pad = 5, top=True, right=True, length=10, width = 1.5)
    plt.title('XZ_Concentrated', fontsize=24)
    plt.rcParams['axes.titley'] = 0.9  
    
for i in range(ct):
    holder[i] = math.radians(holder[i])
    holder_1[i] = math.radians(holder_1[i])
    # hold_q[i] =10*hold_q[i]

lin = lmfit.models.LinearModel()
out1 = lin.fit(holder_1[1:13], x=1/hold_q[1:13])
out2 = lin.fit(holder[1:13], x=1/hold_q[1:13])
# out_holder = np.zeros(len(out2.best_fit))
Lf = 1/out2.best_values['slope']
Lf_1 = 1/out1.best_values['slope']
cept = math.degrees(out2.best_values['intercept'])
cept_1 = math.degrees(out1.best_values['intercept'])

fig, ax = plt.subplots()
plt.scatter(1/hold_q, holder,s=150)
plt.scatter(1/hold_q, holder_1, s=150)
plt.plot(1/hold_q[1:13], out2.best_fit) #out2.best_fit
plt.plot(1/hold_q[1:13], out1.best_fit)
plt.xlabel('1/q (scattering vector, $\mathregular {nm^{-1}}$)', labelpad=20, fontsize=32)
# plt.legend(fontsize=22, frameon=False, loc='best', bbox_to_anchor=(0.6, 0.4))
plt.ylabel('$\mathregular {B_{obs}}$(radians)', labelpad=20, fontsize = 32)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,4])
plt.ylim([0, 10])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.grid(which='major', color='#DDDDDD', linewidth=1.25)
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=1.0)
plt.minorticks_on()
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)
plt.show()

# x = plt.imshow(proc_img_x, cmap = 'gnuplot2', vmin=0, vmax = 1000)
Ruland = pd.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_Compression.xlsx", sheet_name = 'Ruland parameters')
plt.scatter(Ruland.Density, Ruland.Bphi, s = 200)
plt.xlabel('Density (g.$\mathregular {cm^{-3}}$)', labelpad=10, fontsize=30)
# plt.legend(fontsize=22, frameon=False, loc='best', bbox_to_anchor=(0.6, 0.4))
plt.ylabel('Misorientation parameter,° $\mathregular {B_{obs, avg}}$', labelpad=10, fontsize = 30)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,0.13])
plt.ylim([0, 300])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.grid(which='major', color='#DDDDDD', linewidth=1.25)
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=1.0)
plt.minorticks_on()
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)
# def lm_func(x, a1,c1,s1):
#     return a1/(s1*np.sqrt(2*np.pi))*np.exp(((-x-c1)**2)/(2*s1**2))

# gau_ret= lm_func(test1,*out.best_values.values()) #functional response based on parameters
# plt.plot(test1, out.init_fit)
# plt.plot(tth[170:290], Icum[170:290,2])
# plt.plot(tth[170:290], lmfit_Int(Icum[170:290,2], tth[170:290])[0])