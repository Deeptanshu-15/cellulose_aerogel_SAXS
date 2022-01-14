
"""
Created on Tue Jul 20 15:31:01 2021

@author: sid
"""
from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas
import os
import numpy as np
import plotly.express as px
from scipy import integrate, signal
import scipy.optimize as op
import sympy as sy
import lmfit

reader = pandas.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_compression.xlsx", sheet_name='SAXS_intensity_angle')

tth = np.zeros(len(reader))
aero = np.zeros((6,len(reader),2))
# #1 = 0%
# #2 = 85%
# #3 = 40%
# #4 = 70%
# #5 = 98%
# #6 = 95%
#7 = 1.5_0
#8 = 1.5_40
#9 = 2.0_0
#10=0.5_AA
plt.rcParams['figure.dpi'] = 600
plt.rcParams["figure.figsize"] = (10,10)
ctr=0 #column counter
cl=0 #Z-X shifter
for i in range(6):
    for j in range(len(reader)):
        cl=i*3
        for k in range(2):
            aero[i][j] = reader.loc[j, reader.columns[cl]]
            cl+=1
    ctr=cl
    ctr+=1
    # cl=ctr
tth=reader.Theta.values
x_st=-180
x_t=[]
y_st=0
y_t=[]
for i in range(0,375,30):
    x_t.append(x_st+i)  #x_ticks counter and location

for i in range(0,50000,5000):
    y_t.append(y_st+i) #FIgure7a,b

# for i in range(0,11000,1000):
#     y_t.append(y_st+i)

# for i in range(0,5000,1000):
#     y_t.append(y_st+i)

fig, ax = plt.subplots()


def area2(a,z):
    b=np.zeros(len(reader))
    h = a-z
    for i in range(len(reader)-1):
        b[i]=0.5*(h[i+1]-h[i])
    bb=np.absolute(b)
    bs=bb.sum()
    return bs #calculation of the Degree of Alignment based on area between curve and a vs z; then calculating trapezoid area with previous member


def area3(a): #DPO calculation basically area under peak/area under entire curve
    b=np.zeros(len(reader))
    h = a.quantile([0.5])
    holder=[]
    for i in range(len(reader)):
            if (a[i]>h).any():
                holder.append(a[i])
    c=np.zeros(len(holder))
    for i in range(len(reader)-1):
        b[i]=0.5*(a[i+1]-a[i])
    for i in range(len(holder)-1):
        c[i]=0.5*(holder[i+1]-holder[i])
    bb=np.absolute(b)
    cc=np.absolute(c)
    bs=bb.sum()
    cs=cc.sum()
    return cs/bs

def area4(a): # REAL _ DPO calculation basically area under peak/area under entire curve
    b=np.zeros(len(reader))
    h = a.min()
    holder=[]
    denom = a.sum()
    numer = h*360
    return (denom-numer)/denom

# def hermann(a, tth):  #0-360 degree integral
#     num = np.zeros(360)
#     dem = np.zeros(360)
#     for i in range(360):
#             num[i] = a[i]*np.cos(np.deg2rad(tth[i]))**2*np.sin(np.deg2rad(tth[i]))
#             dem[i] = a[i]*np.sin(np.deg2rad(tth[i]))
#     cos = num.sum()/dem.sum()
#     num2=np.absolute(num)
#     dem2=np.absolute(dem)
#     cos2 = num2.sum()/dem2.sum()
#     gamma = 1-(2*cos2)
#     f=(3*gamma-1)/2
#     return f
def hermann(a, tth):
    num = np.zeros(360)
    dem = np.zeros(360)
    for i in range(360):
            num[i] = a[i]*np.cos(np.deg2rad(tth[i]))**2*np.sin(np.deg2rad(tth[i]))
            dem[i] = a[i]*np.sin(np.deg2rad(tth[i]))
    cos = num.sum()/dem.sum()
    num2=np.absolute(num)
    dem2=np.absolute(dem)
    cos2 = num2.sum()/dem2.sum()
    gamma = 1-(2*cos2)
    f=(3*gamma-1)/2
    return f
    
def herman(x,c):
   a=x.max()
   return 4.4 + a/(1+((x-90)/c)**2) + a/(1+((x+90)/c)**2) #lorentzian


#Vivek help for x-axis samples ###WHAT NONSENSE!!!!
# df_reader = pandas.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_compression.xlsx", sheet_name='SAXS_intensity_angle')
# df_a10x_theta = df_reader[['Aerogel_10_x', 'Theta']]
# # df_a10x_theta.plot(x='Theta', y= 'Aerogel_10_x')     
# df_a2z_theta = df_reader[['Aerogel_2_z', 'Theta']]
# # df_a2z_theta.plot(x='Theta', y= 'Aerogel_2_z')
# fift_th_ord_func = np.polyfit(df_a10x_theta['Theta'], df_a10x_theta['Aerogel_10_x'], 6)
# scnd_ord_func = np.polyfit(df_a10x_theta['Theta'], df_a10x_theta['Aerogel_10_x'], 2)
# print("fift_th_ord_func coef : ",fift_th_ord_func )
# print("scnd_ord_func coef : ",scnd_ord_func )
# # fift_th_ord_func_a2z = np.polyfit(df_a2z_theta['Theta'], df_a2z_theta['Aerogel_2_z'], 5)
# # plt.plot(df_a10x_theta['Theta'],df_a10x_theta['Aerogel_10_x'],'o')
# fifth_trendpoly = np.poly1d(fift_th_ord_func) 
# # plt.plot(df_a10x_theta['Theta'],fifth_trendpoly(df_a10x_theta['Theta']), color='#800080', ls = '-.')
# yhat = fifth_trendpoly(df_a10x_theta)
# ybar = sum(df_a10x_theta['Aerogel_10_x'])/len(df_a10x_theta['Aerogel_10_x'])
# SST = sum((df_a10x_theta['Aerogel_10_x'] - ybar)**2)
# SSreg = sum((yhat - ybar)**2)
# R2 = SSreg/SST
# eq = fifth_trendpoly*np.cos(df_a10x_theta['Theta'])**2*np.sin(df_a10x_theta['Theta'])
# res = integrate.quad(eq, 0, np.pi/2)
# eq_2 = fifth_trendpoly*np.sin(df_a10x_theta['Theta'])
# res_2 = integrate.quad(eq_2, 0, np.pi/2)
# cos = res[0]/res_2[0]
# #polynomial regression

# #unique function_scipy curve_fit
# aa, bb = op.curve_fit(lambda t, a, c: a + t*np.sin(t/45+c), tth, df_a10x_theta['Aerogel_10_x'])
# def func(x,a,c):
#     return a + x*np.sin(x/45+c)
# # plt.figure()
# # plt.plot(tth, reader.Aerogel_10_x, 'ko', label="Original Noised Data")
# # plt.plot(tth, func(tth, *aa), 'r-', label="Fitted Curve")
# # plt.legend()
# # plt.show()
# store=np.array(func(tth, *aa))
# eq_a = store*np.cos(np.deg2rad(tth))**2*np.sin(np.deg2rad(tth))
# eq_b = store*np.sin(np.deg2rad(tth))
# cos = eq_a.sum()/eq_b.sum()
# ee=np.absolute(eq_a)
# ee2=np.absolute(eq_b)
# ff2=(3*(ee.sum()/ee2.sum()))/2
# s=s2=0
# for i in range(180):
#     s = s+eq_a[90+i]
#     s2=s2+eq_b[90+i]
# ff = (3*cos-1)/2
# # ff2= (3*(s/s2)-1)/2
# print (ff, ff2)
# # res_a = integrate.quad(eq_a, 0 ,np.pi/2)

# #LMFITparameters
# exp_mod = lmfit.models.ExponentialModel(prefix='exp_')
# pars = exp_mod.guess(reader.Aerogel_1_x, x=tth)

gauss1 = lmfit.models.GaussianModel(prefix='g1_')
pars=gauss1.guess(reader.Aerogel_5_z, x=tth)
# pars.update(gauss1.make_params())

pars['g1_center'].set(value=-150, min=-180, max=-135)
pars['g1_sigma'].set(value=50, min=1)
pars['g1_amplitude'].set(value=100000, min=1000)

gauss2 = lmfit.models.GaussianModel(prefix='g2_')
pars.update(gauss2.make_params())

pars['g2_center'].set(value=10, min=-20, max=40)
pars['g2_sigma'].set(value=50, min=1)
pars['g2_amplitude'].set(value=100000, min=1000)

# mod = gauss1 + gauss2 + exp_mod
mod = gauss1 + gauss2 

init = mod.eval(pars, x=tth)
out = mod.fit(reader.Aerogel_5_z, pars, x=tth)
init = mod.eval(pars, x=tth)
out = mod.fit(reader.Aerogel_5_z, pars, x=tth)

print(out.fit_report(min_correl=0.1))

# fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
# axes[0].plot(tth, reader.Aerogel_1_x, 'b')
# axes[0].plot(tth, init, 'k--', label='initial fit')
# axes[0].plot(tth, out.best_fit, 'r-', label='best fit')
# axes[0].legend(loc='best')
# comps = out.eval_components(x=tth)
# axes[1].plot(tth, reader.Aerogel_1_x, 'b')
# axes[1].plot(tth, comps['g1_'], 'g--', label='Gaussian component 1')
# axes[1].plot(tth, comps['g2_'], 'm--', label='Gaussian component 2')
# # axes[1].plot(tth, comps['exp_'], 'k--', label='Exponential component')
# axes[1].legend(loc='best')

def lm_func(x, a1,c1,s1,a2,c2,s2):
    return a1/(s1*np.sqrt(2*np.pi))*np.exp(((-x-c1)**2)/(2*s1**2)) + a2/(s2*np.sqrt(2*np.pi))*np.exp(((-x-c2)**2)/(2*s2**2))

gau_ret= lm_func(tth,*out.best_values.values()) #functional response based on parameters
holder = holder2 = np.zeros(90)
ss=[]
for i in range(90):
    ss.append(gau_ret[i+180])

# for i in range(90):
#     holder[i] = ss[i]*np.cos(np.deg2rad(i))**2*np.sin(np.deg2rad(i))
#     holder2[i] = ss[i]*np.sin(np.deg2rad(i))

# #curve fitting over
# def hermann(a, tth):
#     num = np.zeros(360)
#     dem = np.zeros(360)
#     for i in range(360):
#             num[i] = a[i]*np.cos(np.deg2rad(tth[i]))**2*np.sin(np.deg2rad(tth[i]))
#             dem[i] = a[i]*np.sin(np.deg2rad(tth[i]))
#     cos = num.sum()/dem.sum()
#     num2=np.absolute(num)
#     dem2=np.absolute(dem)
#     cos2 = num2.sum()/dem2.sum()
#     gamma = 1-(2*cos2)
#     f=(3*gamma-1)/2
#     return f


#PLotting phase


# plt.plot(tth, reader.Aerogel_6_x, color='#FF7F50', ls='--', label='0.102 g.$\mathregular {cm^{-3}}$', dashes =[6,5], linewidth=2.0) #Figure7_c

# plt.plot(reader.Density, reader.Calculated_Area)
# plt.plot(tth, reader.Aerogel_6_z, color='#FF7F50') #Figure7_c

# plt.plot(tth, reader.Aerogel_1_z, label='0.5 wt.% = 9 mg.$\mathregular {cm^{-3}}$', color='#00008B', linewidth=2.0) #Figure7_c
plt.plot(tth, reader.Aerogel_1_x, label='9 mg.$\mathregular {cm^{-3}}$', color='#00008B', linewidth=2.0, dashes = [6,3]) #Figure7_d
# plt.plot(tth, reader.Aerogel_7_z, color='#808080', label='1.5wt.% = 17 mg.$\mathregular {cm^{-3}}$' ,linewidth=2.0)#figure 7_c
# plt.plot(tth, reader.Aerogel_7_x, color='#808080', label='1.5wt.% = 17 mg.$\mathregular {cm^{-3}}$', linewidth=2.0)#figure 7_d
# plt.plot(tth, reader.Aerogel_9_z, color='#753909', label='2.0wt.% = 25 mg.$\mathregular {cm^{-3}}$', linewidth=2.0)#figure 7_c
# plt.plot(tth, reader.Aerogel_9_x, color='#753909', label='2.0wt.% = 25 mg.$\mathregular {cm^{-3}}$', linewidth=2.0) #figure 7_d
plt.plot(tth, reader.Aerogel_3_x, label='46 mg.$\mathregular {cm^{-3}}$', color='#008B8B', dashes =[6,3], linewidth=2.0, ls='--',) #Figure7_a
plt.plot(tth, reader.Aerogel_4_x, label='67 mg.$\mathregular {cm^{-3}}$', color='#800080', dashes =[6,3], linewidth=2.0 , ls='--') #Figure7_a
plt.plot(tth, reader.Aerogel_2_x, label='82 mg.$\mathregular {cm^{-3}}$', color='#FF69B4', dashes =[6,3], linewidth=2.0, ls='--') #Figure7_a
plt.plot(tth, reader.Aerogel_6_x, color='#FF7F50', ls='--', label='102 mg.$\mathregular {cm^{-3}}$', dashes =[6,3], linewidth=2.0) #figure_7_a
# plt.plot(tth, reader.Aerogel_6_x, color='#FF7F50', label='0.102 g.$\mathregular {cm^{-3}}$', dashes =[6,3], linewidth=2.0, ls='--') #figure7_c
plt.plot(tth, reader.Aerogel_5_x, label='112 mg.$\mathregular {cm^{-3}}$', color='#F6A600', dashes =[6,3], linewidth=2.0, ls='--') #Figure7_a
# plt.plot(reader.Aerogel_2_z, tth, label = 'new')
# plt.plot(tth, reader.Aerogel_8_z, label='0.025 g.$\mathregular {cm^{-3}}$', color='#000000', linewidth=2.0, ls='-') #supplement
# plt.plot(tth, reader.Aerogel_8_x, label='0.025 g.$\mathregular {cm^{-3}}$', color='#000000', linewidth=2.0, ls='--') #supplement
# plt.rc('axes', linewidth=2)
# plt.legend(fontsize=22, frameon=False, loc=(0.29,0.52)) #Figure_7a
plt.legend(fontsize=24, frameon=False, loc=(0.6,0.4)) #Figure_7b,c,d
txt = '''XY-Concentrated''' #Header for chart
# plt.xlim([0, 0.15])
# plt.ylim([0, 65000]
plt.xlim([-180,180])
plt.ylim([0,45000])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.xticks(ticks=x_t, labels=x_t, fontsize=24, rotation = 45)
plt.yticks(ticks=y_t, labels=y_t, fontsize=24, rotation = 45)
# plt.ticklabel_format(fontsize = 18)
# plt.ylabel('Intensity (a.U.)', labelpad=20, fontsize=40)
# plt.xlabel('Azimuthal angle, φ (°)', labelpad=20, fontsize = 40)
ax.xaxis.set_tick_params(width=2, size=18)
ax.yaxis.set_tick_params(width=2, size=18)
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)
fig.text(0.6,0.83, txt, fontsize = 26)
plt.show()

new = pandas.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_compression.xlsx", sheet_name='BET data')
new = new.fillna(0)
new = new.replace('-', 0)
plt.style.use('seaborn-white')
x=y=z=np.zeros(len(new))
pore = new.Dpore.tolist()
lamb = new.Thermal_Conductivity.tolist()
strain = new.Strain.tolist()
Int = new.Intensity_z.tolist()    
    
plt.scatter(pore[1:11], lamb[1:11], marker='o', color = 'r', s=150, label = 'Concentration')
plt.scatter(pore[13:19], lamb[13:19], marker = 's', color = '#6495ED', s=150, linewidths=3, label= 'Compression')
plt.scatter(pore[19:22], lamb[19:22], marker = 's', color = '#191970', s=150, linewidths=3)
plt.xlabel('Pore size (nm)', labelpad=20, fontsize=32)
plt.legend(fontsize=22, frameon=False, loc='best', bbox_to_anchor=(0.6, 0.4))
plt.ylabel('Thermal conductivity (mW.$\mathregular {m^{-1}}$.$\mathregular {K^{-1}}$)', labelpad=20, fontsize = 32)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,800])
plt.ylim([15,40])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.grid(which='major', color='#DDDDDD', linewidth=1.25)
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=1.0)
plt.minorticks_on()
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)
plt.show()

dpo = pandas.read_excel(r"C:\Users\sid\Desktop\Cellulose aerogel paper\Lambda_compression.xlsx", sheet_name='DPO')
dpo = dpo.fillna(0)
dpo = dpo.replace('-', 0)
dp = dpo.Actual_DPO.tolist()
her = dpo.Hermann.tolist()
lam = dpo.conductivity.tolist()
dp[9] = 0.13

plt.rcParams["figure.figsize"] = (5,10)
plt.scatter(dp[3:10], lam[3:10], marker='o', color = 'r', s=150, label = 'DPO')
plt.scatter(her[3:10], lam[3:10], marker = 'x', color='b', s=150, label = 'Hermann')
z = np.polynomial.polynomial.polyfit(x, y, 3)
p = np.poly1d(z)
plt.plot(x,p(x),"r--")
plt.xlabel('Fiber orientation (-)', labelpad=20, fontsize=32)
plt.ylabel('Thermal conductivity (mW.$\mathregular {m^{-1}}$.$\mathregular {K^{-1}}$)', labelpad=20, fontsize = 32)
plt.legend(fontsize=22, frameon=False, loc='upper right')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,0.5])
plt.ylim([15,40])
plt.tick_params(axis='both', pad = 10, top=True, right=True)
plt.grid(which='major', color='#DDDDDD', linewidth=1.25)
plt.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=1.0)
plt.minorticks_on()
ax.xaxis.set_tick_params(width=2, size=10)
ax.yaxis.set_tick_params(width=2, size=10)
for axis in ['top', 'bottom', 'left','right']:
    ax.spines[axis].set_linewidth(2.0)
plt.show()