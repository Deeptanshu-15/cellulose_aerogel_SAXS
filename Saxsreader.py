# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 10:03:29 2020

@author: sid
"""

from matplotlib import pyplot as plt
import numpy as np
import plotly.express as px
from scipy import ndimage, signal
import pandas as pd
plt.ion()
#Detector calibration and setup
import pyFAI, pyFAI.detectors, fabio
import pyFAI.distortion as dis
from pyFAI.gui import jupyter
from pyFAI.calibrant import get_calibrant
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import time

#openfile for data #Only RAW DATA has radial averaging, as urwarped is already averaged
fig_z = fabio.open(r"C:\PhD work\PhD_May20\SAXS107cm\Box_01\Aerogel1_2_60s_107cm_01_unwarped.gfrm") #data has additional commands like shape
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
#Plot 1D azimuthal and radial integration--was not removed before; but results are varying. better to use for 1d integrate. radial integration is questionable
#XZ plane plot
res_z = ai.integrate1d(img_z, 2000, unit="q_nm^-1",filename= "integrated.dat", radial_range = [0.3,2.5], correctSolidAngle=False, mask = msk)
rad_z = ai.integrate_radial(img_z, 2000, radial_range= [0,2.5],unit = "chi_deg", correctSolidAngle=True, mask = msk, npt_rad = 200)
jupyter.plot1d(res_z, label = "Compression axis")
jupyter.plot1d(rad_z)

ai.setSPD(SampleDistance = 1.070500031, Center_1=1024.6-.1, Center_2 = 1026.5+0.1) #maybe removed
# ai.setSPD(SampleDistance = 1.070500031, Center_1=xx_ctr, Center_2 = yy_ctr) #weighted average xx_ctr, yy_ctr better for some check and move
d_z = ai.integrate2d(img_z, Desc, correctSolidAngle = False, radial_range = [0.2,2.5], mask = msk, )
# jupyter.plot2d(d_z)

#XY plane plot
res_x = ai.integrate1d(img_x, 2000, unit="q_nm^-1", filename= "integrated.dat", radial_range = [0.3,5], correctSolidAngle=False, mask= msk)
rad_x = ai.integrate_radial(img_x, 500, unit = "chi_deg", radial_range= [0.2,2.5], correctSolidAngle=True, mask=msk)
jupyter.plot1d(rad_x, label = "XY axis")

# #plotting rad same way
# ai.setSPD(SampleDistance = 1.070500031, Center_1=1024.6, Center_2 = 1026.5) #different from above
d_x = ai.integrate2d(img_x, Desc, correctSolidAngle= False, radial_range = [0.2,2.5], mask = msk)
# jupyter.plot2d(d_x)

#Flats plot
# res_flat = ai.integrate1d(flat.data,2000,unit="q_nm^-1",filename="integrated.dat", radial_range = [0.01,10])
# rad_flat = ai.integrate_radial(flat.data,500, radial_range= [0,10],unit = "chi_deg")
# #raw
# res_flat_raw = ai.integrate1d(flat_xz.data,2000,unit="q_nm^-1",filename="integrated.dat", radial_range = [0.01,10])
# rad_flat_raw = ai.integrate_radial(flat_xz.data,500, radial_range= [0,10],unit = "chi_deg", )

#2D Scatter plot without imshow
jupyter.display(img_z)#2d scatter display
jupyter.display(img_x)

#plotting 2d integrated data
intensity, q, tth = d_z
intensity_x, q_x, tth_x = d_x
z_max = np.zeros(len(intensity))
x_max = np.zeros(len(intensity_x))
z_avg = np.zeros(len(intensity))
x_avg = np.zeros(len(intensity_x))

intensity_mask_x = np.zeros(360)
intensity_mask_z = np.zeros(360)


for i in range(len(intensity)): ## Z-axis aligned sample distortion correction for initial beamstop
    z_max[i] = intensity[i].max()
    z_avg[i] = intensity[i].mean()
    x_max[i] = intensity_x[i].max()
    x_avg[i] = intensity_x[i].mean()
    for j in range(Desc):
        if(int(intensity[i][j]>160)):
            intensity_mask_z[i] = j
            break

for i in range(len(intensity_x)): ## X-axis aligned sample distortion correction for initial beamstop

    for j in range(1500):
        if(int(intensity_x[i][j]>160)):
            intensity_mask_x[i] = j
            break
        
#Deletion of Z-axis limits upto mask
corr_intensity_z = np.zeros((360,Desc)) #corrected intensity
corr_intensity_x = np.zeros((360,Desc)) #corrected intensity

for i in range(len(intensity)): #rearrangement for masking Z-axis
    aa = int(intensity_mask_z[i])
    for j in range(Desc):
        if(j>=int(intensity_mask_z[i])):
            k=j-aa
            corr_intensity_z[i][k] = intensity[i][j]
            
for i in range(len(intensity)):
    aa = int(intensity_mask_x[i])
    for j in range(Desc):
        if(j>=int(intensity_mask_x[i])):
            k=j-aa
            corr_intensity_x[i][k] = intensity_x[i][j]

corr_max_x =  np.zeros(len(intensity_x))
corr_avg_x =  np.zeros(len(intensity_x))
corr_sum_x =  np.zeros(len(intensity_x))
corr_max_z =  np.zeros(len(intensity))
corr_avg_z =  np.zeros(len(intensity))
corr_sum_z =  np.zeros(len(intensity))

for i in range(len(intensity)):
    corr_avg_x[i] = corr_intensity_x[i].mean()
    corr_avg_z[i] = corr_intensity_z[i].mean()
    corr_max_x[i] = corr_intensity_x[i].max()
    corr_max_z[i] = corr_intensity_z[i].max()
    corr_sum_x[i] = corr_intensity_x[i].sum()
    corr_sum_z[i] = corr_intensity_z[i].sum()

z_sum = np.zeros(len(intensity))
x_sum = np.zeros(len(intensity_x))
                 
for i in range(len(intensity)): ## Z-axis aligned sample distortion correction for initial beamstop
    z_sum[i] = intensity[i].sum()
    x_sum[i] = intensity_x[i].sum()
    
#moving averaging for correction    
def moving_average(x): #(x,w)
    mv = np.zeros(len(intensity))
    for i in range(len(mv)):
        mv[i] = (x[i-2] + x[i-1] +x[i])/3
    return mv
    # return np.convolve(x, np.ones(w), 'valid') / w

mv_avg = moving_average(x_sum)
# fig_avg = px.scatter(x= tth, y = mv_avg, height = 1200, width = 1200,labels = 'Moving average')
error = (x_sum.max() + x_sum.min())/2
err= np.zeros(len(intensity))
for i in range(360):
    err[i] = abs(np.average(x_sum)-x_sum[i])
print (err.sum(), np.average(x_sum))

index_z = np.argmin(z_sum)
index_x = np.argmin(x_sum) #argumentation for minimum in array

rotate_img_z = ndimage.rotate(img_z, 45, reshape=False) #rotation to index of minimum theta where we can shift/rotate image
rotate_img_x = ndimage.rotate(img_x, abs(tth[index_x]), reshape=False) #rotation to minimum
proc_img_z = rotate_img_z[800:1248, 800:1248]
proc_img_x = rotate_img_x[800:1248, 800:1248]

# rot_z = integrate2d(rotate_img_z, Desc, radialrange=[0.2, 2.5], mask=msk)

fig_int_azi_z = px.scatter(x = tth, y = z_sum, height = 1200, width = 1200,labels = 'Intensity_z', title="Z-axis orientation")
fig_int_azi_z.update_yaxes(title_font=dict(size=24, family='Courier', color='crimson'), title ='Intensity (a.U.)')
fig_int_azi_z.update_xaxes(title_font=dict(size=24, family='Courier', color='crimson'), title ='Azimuthal angle, 𝛘 (°)')
fig_int_azi_z.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
fig_int_azi_z.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
fig_int_azi_x = px.scatter(x = tth_x, y = x_sum, height = 1200, width = 1200, labels = 'Intensity_x', title="X-axis orientation")
fig_int_azi_x.update_yaxes(title_font=dict(size=24, family='Courier', color='crimson'), title ='Intensity (a.U.)')
fig_int_azi_x.update_xaxes(title_font=dict(size=24, family='Courier', color='crimson'), title ='Azimuthal angle, 𝛘 (°)')
fig_int_azi_x.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
fig_int_azi_x.update_yaxes(showline=True, linewidth=1, linecolor='black', mirror=True) #y=corr_sum_x
fig_int_azi_z.show()
fig_int_azi_x.show()

#showing the 2D scatter plot better than jupyter
plt.figure(figsize=(16,9))
plt.axis("off")
z = plt.imshow(proc_img_z, cmap = 'gnuplot2', vmin=0, vmax = 1000)
plt.colorbar() # Show color bar of above image
# patch = patches.Circle((2048, 2048), radius=500)
# x.set_clip_path(patch)
# plt.imsave('98_z', img_z, cmap='gnuplot', dpi = 1200)
plt.show()

x = plt.imshow(proc_img_x, cmap = 'gnuplot2', vmin=0, vmax = 1000)

a,b = res_z
fig, ax = plt.subplots()
q1 = pd.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_Compression.xlsx", sheet_name="q vs Intensity_SAXS")
plt.loglog(q1.Q_range, q1.Aerogel_5_z, label='95% Compressed_XZ', linewidth = 2.0)
plt.loglog(q1.Q_range, q1.Aerogel_5_x, label='95% Compressed_XY', linewidth = 2.0)
plt.loglog(q1.Q_range, q1.Aerogel_1_z, label='0% Compressed_XZ', linewidth = 2.0)
plt.loglog(q1.Q_range, q1.Aerogel_1_x, label='0% Compressed_XY', linewidth = 2.0)
plt.legend(fontsize=22, frameon=False, loc='lower left')
plt.xlabel('Scattering vector, q ($\mathregular {nm^{-1}}$)', labelpad=20, fontsize=32)
# plt.legend(fontsize=22, frameon=False, loc='best', bbox_to_anchor=(0.6, 0.4))
plt.ylabel('Intensity, a.U.', labelpad=20, fontsize = 32)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0.1,5])
plt.ylim([0, 120])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.grid(which='major', color='#DDDDDD', linewidth=1.25)
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=1.0)
plt.minorticks_on()
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)