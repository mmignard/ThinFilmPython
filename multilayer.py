# -*- coding: utf-8 -*-
'''
Created on Tue Feb  7 12:05:30 2023

@author: Marc Mignard
http://www.sspectra.com/sopra.html
https://refractiveindex.info/

'''
import numpy as np
from scipy.optimize import least_squares

def readMaterialFile(fn,w):
    '''
    Read a text file for n&k and interpolate the data at wavelengths w
    '''
    wnk = np.genfromtxt(fn,delimiter='\t',skip_header=5,unpack=True,dtype=float)
    if (np.mean(wnk[0,:])<100): #units are in microns, so convert to nm
        wnk[0,:] = wnk[0,:]*1000
    nk = np.interp(w,wnk[0,:],wnk[1,:]) - 1j*np.interp(w,wnk[0,:],wnk[2,:])
    return np.expand_dims(nk, axis=1)
#wnk = np.genfromtxt('.//materials//SiO2_PVI.txt',delimiter='\t',skip_header=5,unpack=True,dtype=float)

def cauchy(a,w):
    '''
    calculate refractive index at wavelengths w given Cauchy coefficients a[0] & a[1]
    see https://en.wikipedia.org/wiki/Cauchy%27s_equation
    '''

    nk = a[0]+a[1]/np.power(w,2) - 0j; #Cauchy equation assumes non-dispersive material
    return np.expand_dims(nk, axis=1)


def calcCauchyCoeffs(n,w,initial=[1,0]):
    '''
    calcCauchyCoeffs(n,w) - estimates the Cauchy coefficients given the
    measured refractive indices at wavelengths w
    The Cauchy equation works reasonably well for visible wavelengths, but
    the Sellmeier equation is more accurate
    
    Usage: a  = CalcCauchyCoeffs(n,w)

    n = refractive indices     -- dimension (wavelengths X 1 )
    w = free space wavelengths -- dimension (wavelengths X 1 )

    '''
    
    ErrorFunc=lambda a,w,fit: np.squeeze(cauchy(a,w)).real-fit    
    res_lsq = least_squares(ErrorFunc, initial, args=(w, n))

    return res_lsq.x


def sellmeier(a,w):
    '''
    sellmeier(a,w) - returns a vector of refractive indices at the
    
    wavelengths w, given a, the Sellmeier coefficients.
    Usage: n = sellmeier(a,w)
    
    a = Sellmeier coefficients -- dimension (6 X 1 )
    w = free space wavelengths -- dimension (wavelengths X 1 )
    a[4:6] must be in the same units as w (most industry values are in um^2,
    but I usually use nm^2).
    see https://en.wikipedia.org/wiki/Sellmeier_equation
    '''
    wsq = w*w
    n = np.sqrt(1.0 + a[0]*wsq/(wsq-a[3]) + a[1]*wsq/(wsq-a[4]) + a[2]*wsq/(wsq-a[5]))
    return n


def  calcSellmeierCoeffs(n,w,initial = [1,1,1,0,0,0]):
    '''
    calcSellmeierCoeffs(n,w) - estimates the Sellmeier coefficients given the
    measured refractive indices at wavelengths w

    Usage: a  = CalcSellmeierCoeffs(n,w)

    n = refractive indices     -- dimension (wavelengths X 1 )
    w = free space wavelengths -- dimension (wavelengths X 1 )
    
    '''

    ErrorFunc = lambda a,w,fit: np.squeeze(sellmeier(a,w)).real-fit    
    res_lsq = least_squares(ErrorFunc, initial, args=(w, n))

    return res_lsq.x

    
def sqrte(z):
    '''
    evanescent SQRT for waves problems
    Usage: y = sqrte(z)

    z = array of complex numbers
    y = square root of z

    Notes: for z = a-j*b, y is defined as follows:
                [ sqrt(a-j*b),  if b~=0
            y = [ sqrt(a),      if b==0 and a>=0
                [ -j*sqrt(|a|), if b==0 and a<0   (i.e., the negative of what the ordinary SQRT gives)

            this definition is necessary to produce exponentially-decaying evanescent waves
            (under the convention exp(j*omega*t) for harmonic time dependence)

            it is equivalent to the operation y = conj(sqrt(conj(a-j*b))),
            but it fixes a bug in the ordinary SQRT in MATLAB arising whenever the real part is negative
            and the imaginary part is an array with some zero elements. For example, compare the outputs:

            np.conj(np.sqrt(np.conj(-1-0j))) = 0 + 1.000j,
            sqrte(np.array([-1-0j])) = 0 - 1.000j

            np.conj(np.sqrt(np.conj(-1-1j))) = 0.4551 - 1.0987j,
            sqrte(np.array([-1-1j])) = 0.4551 - 1.0987j

    Sophocles J. Orfanidis - 1999-2008 - www.ece.rutgers.edu/~orfanidi/ewa
    Marc Mignard converted from Matlab to Python/numpy
    '''
    y = np.sqrt(z);
    idx = np.where((z.imag==0) & (z.real<0));
    y[idx] = -1j * np.sqrt(np.abs(z[idx]));
    return y

def multilayerR(n,d,w,theta=0,pol='AV'):
    '''
    Usage: refl = multilayerR(n,d,w)
    reflectance of multilayer stack
    adapted from book and Matlab code of Sophocles J. Orfanidis - http://www.ece.rutgers.edu/~orfanidi/ewa

    Parameters
    -------
    n = complex refractive indices -- dimension (wavelengths X materials  )
    d = thicknesses of layers      -- dimension (1           X materials-2)
    w = free-space wavelengths     -- dimension (wavelengths X 1          )
    theta = angle from normal (radians) -- float
    pol = polarization, can be one of
        "TE" = traverse electric or s-polarization
        "TM" = traverse magnetic or p-polarization
        "AV" = average

    Returns
    -------
    refl = reflectance (output)    -- dimension (wavelengths X 1          )
    '''
    Nwavelen = n.shape[0] # number of wavelengths
    M =  n.shape[1]-2     # number of slabs (minus start and end media)
    costh = sqrte(1 - np.power(np.tile(np.reshape(n[:,0]*np.sin(theta),(n.shape[0],1)),(1,M+2))/n,2))

    if ((pol=='TE') or (pol=='AV')):
        nT = n*costh
        r = -np.diff(nT, axis=1)/(nT[:,0:M+1]+nT[:,1:M+2])      # r(i) = (n(i-1)-n(i)) / (n(i-1)+n(i))
        L = np.tile(d,(Nwavelen,1))*n[:,1:M+1]*costh[:,1:M+1]	# n(i)*l(i)*cos(th(i))
        gammaTE = r[:,M]                                        # initialize gamma at right-most interface
        for q in range(M-1,-1,-1):
            delta = 2*np.pi*L[:,q]/w                             #phase thickness in i-th layer
            z = np.exp(-2j*delta)
            gammaTE = (r[:,q] + gammaTE*z) / (1 + r[:,q]*gammaTE*z)

    if ((pol=='TM') or (pol=='AV')):
        nT = n/costh
        r = -np.diff(nT, axis=1)/(nT[:,0:M+1]+nT[:,1:M+2])      # r(i) = (n(i-1)-n(i)) / (n(i-1)+n(i))
        L = np.tile(d,(Nwavelen,1))*n[:,1:M+1]*costh[:,1:M+1]	# n(i)*l(i)*cos(th(i))
        gammaTM = r[:,M]                                        # initialize gamma at right-most interface
        for q in range(M-1,-1,-1):
            delta = 2*np.pi*L[:,q]/w                            #phase thickness in i-th layer
            z = np.exp(-2j*delta)
            gammaTM = (r[:,q] + gammaTM*z) / (1 + r[:,q]*gammaTM*z)

    if (pol=='AV'):
        refl = (np.power(np.abs(gammaTE),2)+np.power(np.abs(gammaTM),2))/2
    else:
        if (pol=='TE'):
            refl = np.power(np.abs(gammaTE),2)
        elif (pol=='TM'):
            refl = np.power(np.abs(gammaTM),2)
        else:
            refl = np.zeros(Nwavelen)
            print('unknown value for "pol"')
    return refl

def multilayerT(n,d,w,theta=0,pol='AV'):
    '''
    Usage: trans = multilayerT(n,d,w)
    transmittance of a multilayer stack
    adapted from book and Matlab code of Sophocles J. Orfanidis - http://www.ece.rutgers.edu/~orfanidi/ewa

    Parameters
    -------
    n = complex refractive indices -- dimension (wavelengths X materials  )
    d = thicknesses of layers      -- dimension (1           X materials-2)
    w = free-space wavelengths     -- dimension (wavelengths X 1          )
    theta = angle from normal (radians) -- float
    pol = polarization, can be one of
        "TE" = traverse electric or s-polarization
        "TM" = traverse magnetic or p-polarization
        "AV" = average

    Returns
    -------
    trans = transmittance (output)  -- dimension (wavelengths X 1          )
    '''
    Nwavelen = n.shape[0] # number of wavelengths
    M =  n.shape[1]-2     # number of slabs (minus start and end media)
    costh = sqrte(1 - np.power(np.tile(np.reshape(n[:,0]*np.sin(theta),(n.shape[0],1)),(1,M+2))/n,2))

    if ((pol=='TE') or (pol=='AV')):
        nT = n*costh
        r = -np.diff(nT, axis=1)/(nT[:,0:M+1]+nT[:,1:M+2])      # r(i) = (n(i-1)-n(i)) / (n(i-1)+n(i))
        t = 2*nT[:,0:M+1]/(nT[:,0:M+1]+nT[:,1:M+2])             # t(i) = 2*n(i-1)/(n(i-1)+n(i))
        L = np.tile(d,(Nwavelen,1))*n[:,1:M+1]*costh[:,1:M+1]	# n(i)*l(i)*cos(th(i))
        tauTE = t[:,M]                                          # initialize tau at right-most interface
        gammaTE = r[:,M]                                        # initialize gamma at right-most interface
        for q in range(M-1,-1,-1):
            delta = 2*np.pi*L[:,q]/w                            #phase thickness in i-th layer
            z1 = np.exp(-1j*delta)
            z2 = np.power(z1,2)
            tauTE = t[:,q]*tauTE*z1 / (1 + r[:,q]*gammaTE*z2)
            gammaTE = (r[:,q] + gammaTE*z2) / (1 + r[:,q]*gammaTE*z2)

    if ((pol=='TM') or (pol=='AV')):
        nT = n/costh
        r = -np.diff(nT, axis=1)/(nT[:,0:M+1]+nT[:,1:M+2])      # r(i) = (n(i-1)-n(i)) / (n(i-1)+n(i))
        t = 2*nT[:,0:M+1]/(nT[:,0:M+1]+nT[:,1:M+2])             # t(i) = 2*n(i-1)/(n(i-1)+n(i))
        L = np.tile(d,(Nwavelen,1))*n[:,1:M+1]*costh[:,1:M+1]	# n(i)*l(i)*cos(th(i))
        tauTM = t[:,M]                                          # initialize tau at right-most interface
        gammaTM = r[:,M]                                        # initialize gamma at right-most interface
        for q in range(M-1,-1,-1):
            delta = 2*np.pi*L[:,q]/w                            #phase thickness in i-th layer
            z1 = np.exp(-1j*delta)
            z2 = np.power(z1,2)
            tauTM = t[:,q]*tauTM*z1 / (1 + r[:,q]*gammaTM*z2)
            gammaTM = (r[:,q] + gammaTM*z2) / (1 + r[:,q]*gammaTM*z2)

    if (pol=='AV'):
        trans = (np.power(np.abs(tauTE),2)+np.power(np.abs(tauTM),2))/2
    else:
        if (pol=='TE'):
            trans = np.power(np.abs(tauTE),2)
        elif (pol=='TM'):
            trans = np.power(np.abs(tauTM),2)
        else:
            trans = np.zeros(Nwavelen)
            print('unknown value for "pol"')
        trans = trans * (n[:,-1]/n[:,0]).real  #eq 5.3.4
    return trans

def xyz2rgb(XYZ):
    '''
    Convert CIE XYZ values to rgb
    '''
    #conversion matrix
    M = np.array([[3.2410,-1.5374,-0.4986],
                [-0.9692,1.8760,0.0416],
                [0.0556,-0.2040,1.0570]])
    rgb = M.dot(XYZ)
    rgb[rgb > 1] = 1
    rgb[rgb < 0] = 0
    
    #find the dark values and save their values
    idx = np.where(rgb<0.00304)
    rgbOld = np.copy(rgb)
    
    #apply gamma conversion
    rgb = np.power(1.055*rgb,1/2.4)-0.055
    
    #ponly do linear conversion for dark values
    rgb[idx] = 12.92*rgbOld[idx]
    
    #clip to range 0..1
    rgb[rgb > 1] = 1
    rgb[rgb < 0] = 0
    return rgb
    
def spec2xyz(spect,illum,cie):
    '''
    Convert a reflecance spectrum to CIE XYZ color coordinates
    spect: array of reflectances at each wavelength
    illum: array of illuminant brightness at each wavelength
    CIE basis vectors, imension 4 x wavelengths (first vector is wavelengths)
    returns XYZ values for each spectrum
    Note: the number of wavelengths in each of spect, illum, and CIE must be the same
    '''
    scaleY = np.sum(illum*cie[2,:])
    XYZ = np.sum(illum*cie[1:4,:]*spect,axis=1)/scaleY
    return XYZ
