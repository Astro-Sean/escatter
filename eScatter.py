#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 09:02:44 2022

@author: seanbrennan



Python implementation of the compton scattering model from 

Comptonization and the shaping of X-ray source spectra - Monte Carlo calculations
-  Pozdnyakov, L. A., Sobol, I. M. & Syunyaev, R. A. 

(https://ui.adsabs.harvard.edu/abs/1983ASPRv...2..189P)

"""


# =============================================================================
# 
# =============================================================================

import time
import numpy as np
import multiprocessing as mp   
import matplotlib.pyplot as plt
from tqdm import tqdm

import os
from numpy.random import choice,shuffle
from itertools import product
import matplotlib.transforms as transforms
from scipy.signal import savgol_filter


# =============================================================================
# 
# =============================================================================
plt.ioff()
# =============================================================================
# Constant
# =============================================================================

c = 299792458 # Speed of light in m s-1
h = 6.626e-34 # Plancks constant in m2 kg s-1
me = 9.10938e-31 # mass of electron in kg
kb = 1.380649e-23 # Stefan boltzmanns constant in m2 kg s-2 K-1
epsilion = 1e-3 #  Probabliliy of photon escape
absorb_probability = 0 # probability that photon is absorbed
sigma_t = 6.652e-29 # thompson cross sectional area m2

# =============================================================================
# Scatterung paramters
# =============================================================================

# Rest wavelength in units of Angstrom
lam_0 = 6563

# Density parameter typoically a number between 0 [Shell like CSM] and 2 [wind like CSM]
s = [2]

# Number of photons to send out
nPhotons = 1e5

# Optical depth (can be a list)
tau_range = [3, 5, 10, 15, 25]

# Electron temperature from just above formation site to edge of boundary
Te_range = [(2e4, 1e4)]

# Inner and outer radius of scattering region (this doesn't change the profile too much)
R_range = [(1e14, 1.01e14)]

# Wind velocity (untested)
vwind = [40]  # km/s

# Shock velocity in km/s
vsh = [3000]


# =============================================================================
# Unit conversion
# =============================================================================

lam_0 = lam_0  * 1e-10 # rest wavelength in meters!!


# =============================================================================
# OUTPUT FOLDERS
# =============================================================================


script_location = os.path.dirname(os.path.abspath(__file__))

# Define folder names
figures_folder = os.path.join(script_location, 'figures')
output_folder = os.path.join(script_location, 'output')

# Create folders if they don't exist
if not os.path.exists(figures_folder):
    os.makedirs(figures_folder)
    print(f"Created 'figures' folder at: {figures_folder}")
else:
    print(f"'figures' folder already exists at: {figures_folder}")

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created 'output' folder at: {output_folder}")
else:
    print(f"'output' folder already exists at: {output_folder}")
    
        
        
# =============================================================================
# Functions for convience
# =============================================================================

lam_2_nu = lambda lam: c/lam #  Wavelength to frequency
nu_2_lam = lambda nu: c/nu # Frequency to wavelngth
E = lambda nu: nu * h # photon energy

# =============================================================================
# Setup
# =============================================================================


def calculate_bins(x):
    # https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
    # https://stackoverflow.com/questions/33203645/how-to-plot-a-histogram-using-matplotlib-in-python-with-a-list-of-data
    import numpy as np
    q25, q75 = np.nanpercentile(x, [25, 75])
    bin_width = 2 * (q75 - q25) * len(x) ** (-1/3)
    bins = round((np.nanmax(x) - np.nanmin(x)) / bin_width)
    return bins



def lambda2velocity(lam,lam0):
    
    "Convert wavelength shift to velocity in km/s"

    v = c*((lam-lam0)/lam0)
    return v/1000

def calculate_A(tau, r1, r2, s):
    """
    Calculate the constant A given the Compton optical depth tau, 
    the Thomson cross-section sigma_T, the distances r1 and r2, and the exponent s.
    
    Parameters:
    tau (float): Compton optical depth
    sigma_T (float): Thomson cross-section
    r1 (float): Initial distance
    r2 (float): Final distance
    s (float): Exponent in the electron number density equation
    
    Returns:
    float: The constant A
    """

    if np.isclose(s,1):
        A = tau / (sigma_t * np.log(r2 / r1))
    else:
        A = (tau * (1 - s)) / (sigma_t * (r2**(1 - s) - r1**(1 - s)))
    
    return A


def velocity2lambda(v,lam0):

       
    "Convert velocity shift to wavelength"
    lam0 =  ((v * lam0 )/c) + lam0
    
    return lam0


def BoundaryDistance(r,omega,R=1):
    """
    Compute distance for origin (r = {0,0,0}) to boundary of sphere at R = 1
    """
    l = -1 * np.dot(r,omega) + (R**2 - np.dot(r,r) + np.dot(r,omega)**2)**0.5
    return l


def EscapeProbability(l,mfp):
    """Probably that photon escape the scattering region"""
    P = np.exp(-l/mfp)
    return P 



def Ne(tau,R1, R2):
    """Number of electrons given a certain optical depth"""
    return tau/(sigma_t * (R2 - R1))


def mu_prime(v_div_c):
    e = np.random.random()
    mu_prime = (v_div_c + (2*e) - 1)/(1 + ((v_div_c)*(2*e - 1)))
    return mu_prime


def phi():
    """Pick a scattering angle following a sin^2(x) distribution"""
    e2 = np.random.random()
    draw =2*np.pi*e2
    return draw

def re():
    e = 1.602e-19  # Elementary charge in Coulombs
    pi = np.pi  # Pi constant
    epsilon_0 = 8.854e-12  # Vacuum permittivity in C^2/(N*m^2)

    
    # Calculate electron radius
    r_e = e**2 / (4 * pi * epsilon_0 * me * c**2)
    return r_e


def rho(v):
    """rho parameter from section 9.5 """
    return np.sqrt(v[0]**2 + v[1]**2)


def omega_prime(mu_prime,v,phi,rho):
    
    
    """
    Direction Vector after the scattering event
    """
    
    omega_1_prime = (mu_prime*v[0]) + np.sqrt(1-mu_prime**2) * (1/rho) * (   v[1]*np.cos(phi) + v[0]*v[2]*np.sin(phi))
    
    omega_2_prime = (mu_prime*v[1]) + np.sqrt(1-mu_prime**2) * (1/rho) * (-1*v[0]*np.cos(phi) + v[1]*v[2]*np.sin(phi))
    
    omega_3_prime = (mu_prime*v[2]) - np.sqrt(1-mu_prime**2) * (rho) * np.sin(phi)

    omega_prime = np.array([omega_1_prime,omega_2_prime,omega_3_prime])
    
    return omega_prime


def x_prime(x,alpha,nu,gamma,mu_prime,v_div_c):
    
    """
    x' paramter from section 9.5, see equation 2.8
    """
    ratio = 1 +  ((h*nu*(1-alpha) )) / ( (gamma*me*(c**2))*(1 - (mu_prime * v_div_c)))

    return x * 1/ratio 


def sigma_hat(x):
    
    """
    Cross sectional approximaton from Section 9.4
    """
    
    if x<=0.5:
        return (1/3) + (0.141*x) - (0.12*x**2) + (1 + 0.5*x)*((1+x)**-2)
        
    elif x>0.5 or x<=3.5:
        return (np.log(1+x) + 0.06) * (1/x)
        
    else:
        return (np.log(1+x) + 0.5 - (2 + 0.076*x)**-1) * (1/x)
    
def X(x,xp):
    return (x/xp) + (xp/x) + 4*(1/x - 1/xp) + 4*(1/x - 1/xp)**2




def sample_annulus(inner_radius = 1,outer_radius=2):


    
    r = np.sqrt(np.random.uniform(inner_radius**2, outer_radius**2))
    theta = np.random.uniform(0, 2*np.pi)
    
    # Convert polar coordinates to Cartesian coordinates
    x =  r*np.cos(theta)
    y =  r*np.sin(theta)
    z =  np.random.uniform(-1, 1)
    
    return np.array([x,y,z])



def sample_spherical(r):

    theta = 2*np.pi*np.random.random()
    phi =     np.pi*np.random.random()
    
    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    
    return np.array([x,y,z])



def sample_hemisphere(r=1):

    theta = np.pi*np.random.random()
    phi =   np.pi*np.random.random()
    
    
    x = r * (2*(np.cos(theta) * np.sin(phi))-1)
    y = r * (2*(np.sin(theta) * np.sin(phi))-1)
    z = r * (2*(np.cos(phi))-1)
    return np.array([x,y,z])



def Te_func(r,T_in,T_out,R1,R2):

    
    # T = T_in*(1/(1+r))**0.4
    
    m = ( T_out - T_in) / (R2 - R1)
    c = T_in - m * R1
    
    T_r = m * r + c
    return  T_r


vrange = np.arange(0,c,1000000)

# =============================================================================
# 
# =============================================================================

def maxwellian_velocity_probability(v, T):

    prefactor = np.sqrt((me / (2 * np.pi * kb * T))**3)
    exponent = - (me * v**2) / (2 * kb * T)
    
    return prefactor * 4 * np.pi * v**2 * np.exp(exponent)




def get_electron_velocity(T):

    prob = maxwellian_velocity_probability(vrange, T)
    prob/=np.sum(prob)
    v = choice(vrange ,1, p = prob)[0]
    
    return v


# =============================================================================
# Mean free Path
# =============================================================================

def MFP(nu,Ne,Te):
    
    
    if h*nu < me*c**2 and kb*Te < me*c**2: 
        l_inv = Ne*sigma_t*(1 - 2 * (h*nu / (me*c**2)) - 5*(h*nu / (me*c**2)) * (kb*Te / (me*c**2)))
        
    elif h*nu*kb*Te < (me*c**2)**2 and kb*Te > me*c**2:
    
        l_inv = Ne*sigma_t*(1 - 8 *(h*nu/(me*c**2))*(kb*Te/(me*c**2)))

    elif h*nu > me*c**2 and kb*Te > me*c**2:
        
        l_inv = 3/16 * (Ne*sigma_t) * ((me*c**2)/(h*nu)) * ((me*c**2)/(kb*Te)) * ( 0.077 + np.log(4 * (h*nu / (me*c**2)) * (kb*Te / (me*c**2))))

    elif h*nu > me*c**2 and kb*Te < me*c**2:
        
        l_inv = 3/8 * (Ne*sigma_t) * ((me*c**2)/(h*nu))* (np.log( (2 * h*nu / (me*c**2)) ) + 0.5 +   (kb*Te / (me*c**2))) * (1-(3/2)*(kb*Te / (me*c**2)))
        
        
    else:

        raise Exception('somthing is broken with the MFP - contact someone!')
        
    return 1/l_inv

 
# =============================================================================
# 
# =============================================================================
import random


def fwhm_kms_to_angstrom(rest_wave, fwhm_kms):
    # Speed of light in km/s
    speed_of_light_kms =c/1000
    
    # Convert FWHM from km/s to Ångström
    fwhm_angstrom = rest_wave * (fwhm_kms / speed_of_light_kms)
    
    return fwhm_angstrom





def gaussian_input(rest_wave,FWHM_kms):
    
    # FWHM_ms = FWHM_kms * 1000
    fwhm_lam = fwhm_kms_to_angstrom(fwhm_kms=FWHM_kms ,rest_wave  = lam_0)
    lam =  random.gauss(rest_wave, sigma = fwhm_lam / 2.235)
    return lam

# =============================================================================
# 
# =============================================================================

def scatterPhoton(args):

    
    photon_absorbed = False
    
    
    # For each photon interacting, grab a random electron
    T_in  = args[0][0]
    T_out  = args[0][1]
    tau = args[1]

    s = args[6]
    # print(s)
    R1   = args[3][0]
    R2   = args[3][1]
    
    # Convert from cm to m
    R1/=100
    R2/=100
    
    vwind = args[4]
    vshock = args[5]
    

    # RANDOM DIRECTION OF PHOTON (isotropic)
    omega0 = sample_spherical(r = 1)
    
    # Position where the photon is created
    segs = np.linspace(R1*1.001,R2,10000)  
    segs_prob = 1/((segs)**4)
    r =   choice(segs ,1,p=segs_prob/np.sum(segs_prob))[0]

    r0 = sample_spherical(r = r)
    
    rMag = np.sqrt(sum([i**2 for i in r0]))
    

    # V =  -1 * 1000 * vshock *((R1/rMag)**2)* (r0/np.linalg.norm(r0)) - vwind*1000
    
    V = -1 * 1000 * vshock * r0/np.linalg.norm(r0)
    lam_random = velocity2lambda(V[0],lam_0)
    # lam_random = gaussian_input(rest_wave = lam_0,FWHM_kms = vshock)
    

    nu0 = lam_2_nu(lam_random)
    

    # Probability of escape
    w0 = 1
    
    # Scattering counter
    scatters = 0
    
    # Normalisation constant for r^-s
    # A = tau/sigma_t  *  ((1/R1) - (1/R2))**-(s-1)
    # print(A)

    # A = tau/sigma_t  * ((-s + 2)/(R2**(-s+2) - R1**(-s+2)))
    A = calculate_A(tau, R1, R2, s)

    
    
    # Denisty as a function of distance from inner boundary
    Ne_density = lambda r: A * r **-s

    
    pathTaken = []
    photonWeights = []
    photonNu = []
    

    while True:
        
        while True:

            rMag = abs(np.sqrt(sum([i**2 for i in r0])))
          
            #  Find a random electron 
            e1 = np.random.random()
            e2 = np.random.random()
            
            #  Get its velocity
            v3_e = 2*e1 - 1
            v2_e = np.sqrt(1 - v3_e**2) * np.sin(2*np.pi*e2) 
            v1_e = np.sqrt(1 - v3_e**2) * np.cos(2*np.pi*e2) 

            v_e = np.array([v1_e,
                            v2_e,
                            v3_e]) 
            rho1 = rho(v_e)
            

            # Temperature of current region
            Te_r = Te_func(rMag,T_in,T_out,R1,R2)
            
            # Find the photons distance to boundary
            l0 = BoundaryDistance(r0,omega0,R2)
            
            # Mean free path for travelling photon
            mfp = MFP(nu0,Ne =  Ne_density(rMag) ,Te = Te_r)
                
            # PROBABILITY THAT A PHOTON WILL ESCAPE
            Li = EscapeProbability(l0, mfp)
           
            # Weighting each photon by it's probability that it escapes
            w1  = w0 - w0*Li
            
            #  Get the MFP
            lambda_i = -1 * mfp * np.log(1 - np.random.random()*(1 - Li))
            

            # Momentum of Electron
            vMag = get_electron_velocity(Te_r)
            p = vMag*me
            eta = p/(me*c)
        
            # Lorentz factor    
            gamma = np.sqrt((eta**2) + 1)

            v_div_c = eta / gamma
            mu0 = np.dot(v_e,omega0)
        
        
            x = 2 * gamma * ((h * nu0) / (me*c**2)) * (1 - mu0*v_div_c )
        
            while True:
                
                #  Get posible direction of the scattering 
                mu1 = mu_prime(v_div_c)
                phi1 = phi()
                
                omega1 = omega_prime(mu1,
                                     v_e,
                                     phi1,
                                     rho1)
                
                #  deflection angle
                alpha = np.dot(omega0,omega1) 
    
                x1 = x_prime(x,alpha,nu0,gamma,mu1,v_div_c)
                
                Y = ((x1/x)**2) * X(x,x1)

                
                if  2 * np.random.random() < Y : 
                    # Found a good scattering event
                    break
        
            # Scattered Wavelength 
            nu1 = (1/h) * (x1) *  (1 / (2*gamma*(1 - (mu1*(v_div_c))))) * me * (c**2)
            
            
            # Accepted the scattering event
            if np.random.random() < 0.375*sigma_hat(x)*(1 - mu0*v_div_c): 
                
                       
                # if absorb_probability > np.random.random() : 
                #     photon_absorbed = True
                    

                #     break
                break
        
 
                
            
        # Follow to photon from later plotting
        pathTaken.append((r0[0],r0[1],r0[2]))
        
        #  New position vector on photon
        r0 = r0 + lambda_i*omega0
        
        #  New position scalar
        rMag = abs(np.sqrt(sum([i**2 for i in r0])))
        
        
        #  If the photon has gone backwards or something funky has happened - ignore the photon
        if np.isnan(rMag) or rMag < R1  :            
            
            photonNu = [np.nan]
            photonWeights = [np.nan]
            scatters = np.nan
            pathTaken = [np.nan]
            
            break
        
               
        weight = w0*Li
        # if photon_absorbed:
        #     omega1 = sample_spherical(r = 1)
        #     nu1 = nu0
        #     weight = 0
        #     photon_absorbed = False
            
            
        # Append the photon details
        photonWeights.append(weight)
        photonNu.append(nu0)
        
        
        # The photon has escaped
        if w1<epsilion: 
            
            # Ignore photons travelling in the negative x direction (occulaion by the photosphere) 
            if ( r0[0] < 0):
                photonNu = [np.nan]
                photonWeights = [np.nan]
                scatters = np.nan
                pathTaken = [np.nan]
                
            break


        

        scatters+=1
   
 



        
        w0 = w1
        omega0 = omega1
        nu0 = nu1
        
        
    return scatters, photonNu,photonWeights,pathTaken
       
    
    

# =============================================================================
#  Run with multiprocessing handler
# =============================================================================


def mp_handler(Te=10000, tau=10, tau_CSM=None, R=100, vwind=200, vshock=2000, 
               nPhotons=1e6, s=2, returnPath=False, randomselection=False, photonHist=None):
    lams = []
    scatters = []
    pathTaken = []
    photonWeights = []

    start_time = time.time()

    nCPU = 8  # Default number of CPUs to use
    chunksize = int(nPhotons) // (4 * nCPU)
    extra = int(nPhotons) % (4 * nCPU)
    print(f'Running model on {nCPU} CPUs\n')

    if extra:
        chunksize += 1

    args = [[Te, tau, tau_CSM, R, vwind, vshock, s, returnPath, randomselection, photonHist]] * int(nPhotons)

    with mp.Pool(processes=nCPU) as pool:
        for res in tqdm(pool.imap_unordered(func=scatterPhoton, iterable=args), total=int(nPhotons)):
            scatters.append(res[0])
            lams.extend(res[1])
            photonWeights.extend(res[2])
            pathTaken.append(res[3])

        pool.close()
        pool.join()

        for p in pool._pool:
            if p.exitcode is None:
                p.terminate()

    print('\n\tTime taken: %.1fs' % (time.time() - start_time))
    return scatters, lams, photonWeights, pathTaken


    
# =============================================================================
# Get parameters ready
# =============================================================================



modelParams = list(product(tau_range,Te_range,R_range,vsh,vwind,s))
shuffle( modelParams )

args = modelParams[0]



print('Total number of models: %d' % len(modelParams))


# =============================================================================
# 
# =============================================================================

plt.ioff()
plt.close('all')

returnPath = False
redo = True

# =============================================================================
# 
# =============================================================================


if __name__ == '__main__':
           
    counter = 0
    fig = plt.figure()
    
    ax1 = fig.add_subplot(121)
    ax2 = ax1.twiny()
    ax3 = fig.add_subplot(122)
    
    for i in modelParams:
        

             
 
        nl = '\n'


        tau = round(i[0],1)
        Te = i[1]
        R = i[2]
        vshock = i[3]
        vwind = i[4]    
        s = i[5]
        start = time.time()
        
        counter+=1
        
        fname = fr'LAM0_{int(lam_0/1e-10)}_tau_{tau}_T_e_{Te}_R_{R}_vs_{vshock}_vw_{vwind}_s_{s}.txt'
        fpath = os.path.join(output_folder, fname)
        
        if os.path.exists(fpath) and not redo:
            print('Skipping %s as file already exists' % fname )
            continue
        
        
        print('\n[%d / %d]' % (counter,len(modelParams)))
        
        mH = 1.67262192e-27      # hydrogen mass in kg
        msun = 1.989e+30 # solar mass in kg
        
        Ne_total =  Ne(tau,R1 = R[0]/100,R2 = R[1]/100)
        mass_material = (mH * Ne_total) * 4/3 * np.pi * ( (R[1]/100)**3 - (R[0]/100)**3)
        
        text = f"""
       Parameters:
            
       Scattering {nPhotons:.0e} photons through a shell of R1 = {R[0]:.1e} cm and R2 = {R[1]:.1e} cm
       with Te = {Te[0]:.0f} K at the inner boundary and {Te[1]:.0f} K at the outer.
        
       The model has τ = {tau}, which gives an approximate mass of {mass_material/msun:.1e} M☉.

        """
        print(text)

        
        results = mp_handler(Te = Te,
                             tau = tau,
                             R = R,
                             nPhotons = nPhotons,
                             vwind = vwind,
                             vshock = vshock,
                             s = s,
                             returnPath  = returnPath )

        scatters,photonNu,photonWeights,pathTaken = results

        
        output_photons = [i for i in pathTaken if isinstance(i[0],tuple)]
        print('\n\t%.1f%% photons made it out alive' % (100*len(output_photons)/nPhotons))
        
        
        
        
        photonNu = np.array(photonNu)
        
        idx = ~np.isnan(photonNu)
        photonWeights = np.array(photonWeights)[idx]
        photonNu = photonNu[idx]
        
        
        # idx = photonWeights>0
        
        # photonWeights  = photonWeights[idx]
        # photonNu = photonNu[idx]

        photonLams = (c/photonNu)/1e-10

        dl = 1
        
        photonLams_bins = np.linspace(min(photonLams),max(photonLams),1000)
        
        
        flx, photonLams_bins = np.histogram(photonLams,bins = photonLams_bins )
        photonLams_bins = np.array([(photonLams_bins[i]+photonLams_bins[i+1])/2. for i in range(len(photonLams_bins)-1)])
        
            
        
        flx_smooth  = savgol_filter(flx, 10,2)
        
        
        norm = np.max(flx_smooth[abs(photonLams_bins- lam_0/1e-10 < 5)])
        
   
        flx_smooth =  flx_smooth / norm
        flx = flx /norm
        
        text = '$\\tau$=%.1f \n$V_{shock}$ = %.0f \n$R_{1/2}$ = %.1f' % (tau,vshock,R[1]/R[0])
        
        
        l = ax1.step(photonLams_bins- lam_0/1e-10,
                  flx_smooth,label = text,lw = 0.5)
        
        ax1.step(photonLams_bins- lam_0/1e-10,
                  flx,lw = 0.5,alpha = 0.05,color = l[0].get_color())
        
        ax1.step(-1*(photonLams_bins- lam_0/1e-10),
                 flx_smooth,color = 'grey',ls = '--',lw = 0.25)


# =============================================================================
# 
# =============================================================================
        

        fpath = os.path.join(output_folder, fname)
              
        np.savetxt(fpath, np.vstack([photonLams_bins,flx/np.nanmax(flx)]).T,
                    delimiter=', ')
    

# =============================================================================
    # plot velocity axis
# =============================================================================


        
        mn, mx = ax1.get_xlim()

           
        ax2.set_xlim(xmin = lambda2velocity(mn+lam_0/1e-10,lam0 = lam_0/1e-10),
                     xmax = lambda2velocity(mx+lam_0/1e-10,lam0 = lam_0/1e-10))
        

        
        ax2.axvline( -1 *vshock,lw = 0.5,ls = '--',color = l[0].get_color(),alpha = 0.25)
  
        
        ax1.axvline(0 ,lw = 0.5,ls = ':',color = 'grey')
        ax1.legend(frameon = False,fontsize = 7)
        
                 
      
# =============================================================================
#      Histogram with number of scatters   
# =============================================================================

        scatters = np.array(scatters)[np.isfinite(scatters)]
        
        counts, bins = np.histogram(scatters )
        
        ax3.hist(scatters,
                 # bins = calculate_bins(scatters),
                  histtype = 'step',lw = 0.5,color = l[0].get_color())
        
        ax3.axvline(tau**2,label ='$\\tau^2$',ls = '--',
                    color = l[0].get_color(),lw = 0.5)
 
        
        trans = transforms.blended_transform_factory(ax3.transData, ax3.transAxes)
        ax3.text(tau**2,0.01,'$\\tau_e^2$',rotation = 90,transform=trans,
                 ha = 'left',va = 'bottom',color = l[0].get_color())
        
        
        handles, labels = ax1.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        
# =============================================================================
#         Save the output profile
# =============================================================================

        
        ax1.set_xlabel('$\lambda~[A]$')
        ax2.set_xlabel('Velocity [$km~s^{-1}$]')
        ax3.set_xlabel('Number of\nScattering events')
        

        ax3.set_yscale('log')

        
        ax1.set_ylabel('$Normalised Flux$')
        ax3.set_ylabel('$Occurance$')
        
        fig.tight_layout()
        
        fpathFig = os.path.join(figures_folder, fr'tau_{tau}_T_e_{Te}_R_{R}_vs_{vshock}_vw_{vwind}_s_{s}.png')
        plt.savefig(fpathFig,dpi = 300)
        
      
                
                
        plt.show()
        
    
# =============================================================================
#     
# =============================================================================

                            
                                 
                    