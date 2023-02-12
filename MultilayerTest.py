"""
Created on Wed Oct 16 09:57:10 2019
Used to verify operation of code in multilayer.py
"""

import matplotlib.pyplot as plt
import numpy as np
import multilayer as ml

saveImages = False
fpath = './/media//'

##########################################################################
###    Fit Cauchy coefficients to glass index of refraction 
###     
##########################################################################
w = np.linspace(360,830,95)
fn = 'SiO2_PVI'
#fn = '7059'
n = np.squeeze(ml.readMaterialFile('.//materials//'+fn+'.txt',w)).real
a = ml.calcCauchyCoeffs(n,w)
b = ml.calcSellmeierCoeffs(n,w)

plt.figure()
plt.subplot()
plt.plot(w,n,label='measured data')
plt.plot(w,ml.cauchy(a,w),label='Cauchy fit data')
plt.plot(w,ml.sellmeier(b,w),label='Sellmeier fit data')
plt.legend()
plt.grid(True)
plt.title('Fit Cauchy coefficients to glass index of refraction\nfile='+fn)
plt.xlabel('wavelength (nm)')
plt.ylabel('index of refraction')
#plt.xlim([theta[0],theta[-1]])
#plt.ylim([0,1])
plt.show

##########################################################################
###    reflection from air/glass interface 
###     
##########################################################################

w = np.array([650])
d = np.array([])    #thickness of layers, first and last layers are infinite
n = np.ones((np.size(w),1)) #air
n = np.concatenate((n,ml.cauchy([1.5220,0.00459],w)),axis=1) #crown glass
theta = np.linspace(0,89.999,50) #div by 0 error at 90
reflTE = np.zeros(np.size(theta))
reflTM = np.zeros(np.size(theta))
for t in np.arange(np.size(theta)):
    reflTE[t] = ml.multilayerR(n,d,w,theta[t]*np.pi/180,'TE')
    reflTM[t] = ml.multilayerR(n,d,w,theta[t]*np.pi/180,'TM')
plt.figure()
plt.subplot()
plt.plot(theta,reflTE,label='refl TE, s-pol')
plt.plot(theta,reflTM,label='refl TM, p-pol')
plt.legend()
plt.grid(True)
plt.title('Reflectance from air to glass at {}nm'.format(w[0]))
plt.xlabel('angle from normal (degrees)')
plt.ylabel('reflectance')
plt.xlim([theta[0],theta[-1]])
plt.ylim([0,1])
if saveImages:
    plt.savefig(fpath+"Air2Glass.svg", bbox_inches='tight')
plt.show

##########################################################################
###    reflection from 1/4 wave anti reflection coating 
###     
##########################################################################

w = np.array([650])
d = np.array([650/1.38/4])    #thickness of layers, first and last layers are infinite
n = np.ones((np.size(w),1)) #air
#n = np.concatenate((n,ml.cauchy([1.23,0.00459],w)),axis=1) #make believe
n = np.concatenate((n,ml.readMaterialFile('.//materials//MgF2_sopra.txt',w)),axis=1)
n = np.concatenate((n,ml.cauchy([1.5220,0.00459],w)),axis=1) #crown glass
theta = np.linspace(0,89.999,50) #div by 0 error at 90
reflTE = np.zeros(np.size(theta))
reflTM = np.zeros(np.size(theta))
for t in np.arange(np.size(theta)):
    reflTE[t] = ml.multilayerR(n,d,w,theta[t]*np.pi/180,'TE')
    reflTM[t] = ml.multilayerR(n,d,w,theta[t]*np.pi/180,'TM')
plt.figure()
plt.subplot()
plt.plot(theta,reflTE,label='refl TE, s-pol')
plt.plot(theta,reflTM,label='refl TM, p-pol')
plt.legend()
plt.grid(True)
plt.title('Reflectance from air to glass at {}nm'.format(w[0]))
plt.xlabel('angle from normal (degrees)')
plt.ylabel('reflectance')
plt.xlim([theta[0],theta[-1]])
plt.ylim([0,1])
plt.show

##########################################################################
###    Mirasol display (interferometric modulator) 
###     
##########################################################################
w = np.linspace(360,830,95)
d = np.array([5,50,250])    #thickness of layers, first and last layers are infinite
n = ml.readMaterialFile('.//materials//SiO2_PVI.txt',w);
n = np.concatenate((n,ml.readMaterialFile('.//materials//MoCr_PVI.txt',w)),axis=1)
n = np.concatenate((n,ml.readMaterialFile('.//materials//SiO2_PVI.txt',w)),axis=1)
n = np.concatenate((n,np.ones((np.size(w),1))),axis=1) #air
n = np.concatenate((n,ml.readMaterialFile('.//materials//Al_PVI.txt',w)),axis=1)

plt.figure()
plt.subplot()
for tAir in [30,150,200,250]:
    d[2] = tAir    
    plt.plot(w,ml.multilayerR(n,d,w),label='tAir = {}nm'.format(tAir))
plt.legend()    
plt.grid(True)
plt.title('Reflectance of variable airgap in interferometric display')
plt.xlabel('wavelength (nm)')
plt.ylabel('reflectance')
plt.xlim([w[0],w[-1]])
plt.ylim([0,1])
if saveImages:
    plt.savefig(fpath+"imodSpectrum.svg", bbox_inches='tight')
plt.show
plt.show

cie = np.genfromtxt('ciexyz31_5.txt',delimiter=',',skip_header=0,unpack=True,dtype=float)
w = cie[0,:]
D = np.genfromtxt('D65spect.txt',delimiter=',',skip_header=0,unpack=True,dtype=float)
D65 = np.interp(w,D[0,:],D[1,:])

inc = 2
dr = np.arange(0,500,inc)
nc = dr.size #number of colors
c = np.arange(0,nc,1)
fig = plt.figure(figsize=(5,2),dpi=300)
ax = fig.add_subplot(111)
for p in c:
    d[2] = dr[p]
    rgb = ml.xyz2rgb(ml.spec2xyz(ml.multilayerR(n,d,w),D65,cie))
    ax.add_patch(plt.Rectangle((inc*p, 0), inc, 1,color=(rgb[0],rgb[1],rgb[2])))
plt.xlim([0,d[-1]])
plt.ylim([0,1])
plt.yticks([])
plt.title('Color of interferometric display element')
plt.xlabel('airgap thickness (nm)')
if saveImages:
    plt.savefig(fpath+"imodColors.svg", bbox_inches='tight')
plt.show
plt.show


##########################################################################
###    colors of titanium dioxide layer on titanium 
###     
##########################################################################

cie = np.genfromtxt('ciexyz31_5.txt',delimiter=',',skip_header=0,unpack=True,dtype=float)
w = cie[0,:]
D = np.genfromtxt('D65spect.txt',delimiter=',',skip_header=0,unpack=True,dtype=float)
D65 = np.interp(w,D[0,:],D[1,:])
d = np.array([5])    #thickness of layers, first and last layers are infinite
n = np.ones((np.size(w),1)) #air
n = np.concatenate((n,ml.readMaterialFile('.//materials//TiO2_sopra.txt',w)),axis=1)
n = np.concatenate((n,ml.readMaterialFile('.//materials//Ti_sopra.txt',w)),axis=1)

inc = 1
d = np.arange(0,180,inc)
nc = d.size #number of colors
c = np.arange(0,nc,1)
fig = plt.figure(figsize=(5,2),dpi=300)
ax = fig.add_subplot(111)
for p in c:
    rgb = ml.xyz2rgb(ml.spec2xyz(ml.multilayerR(n,[d[p]],w),D65,cie))
    ax.add_patch(plt.Rectangle((inc*p, 0), inc, 1,color=(rgb[0],rgb[1],rgb[2])))
plt.xlim([0,d[-1]])
plt.ylim([0,1])
plt.yticks([])
plt.title('Color of thin oxide on titanium')
plt.xlabel('oxide thickness (nm)')
if saveImages:
    plt.savefig(fpath+"TiOx_color.svg", bbox_inches='tight')
plt.show
plt.show

