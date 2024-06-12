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

# =============================================================================
# Constant
# =============================================================================

c = 299792458 # Speed of light in m s-1
h = 6.626e-34 # Plancks constant in m2 kg s-1
me = 9.10938e-31 # mass of electron in kg
kb = 1.380649e-23 # Stefan boltzmanns constant in m2 kg s-2 K-1
epsilion = 1e-1 #  Probabliliy of photon escape
sigma_t = 6.652e-29 # thompson cross sectional area m2

# =============================================================================
# Scatterung paramters
# =============================================================================

# rest wavelngth in units of Angstrom 
lam_0 = 6563

# denisty parameter either 0 or 2 [0 untested]
s = [2]

# Number of photons to send out
nPhotons = 1e4

# Optical depth (can be list)
tau_range = [5]

# Electron temperature from just above formation site to edge of boundary
Te_range = [(30e3,10e3)]

#unitless radius of scattering region
R_range = [10]



# wind velocity (untested)
vwind = [0] #km/s


# Shock velocity in km/s 
vsh = [5500]


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

def lambda2velocity(lam,lam0):
    
    "Convert wavelength shift to velocity in km/s"

    v = c*((lam-lam0)/lam0)
    return v/1000


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



def Ne(tau,R = 1):
    """Number of electrons given a certain optical depth"""
    return tau/(sigma_t * R)


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



def Te_func(r,T_in):

    
    T = T_in*(1/(1+r))**0.4

    return  T


vrange = np.arange(0,c,100000)

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

def scatterPhoton(args):

    # For each photon interacting, grab a random electron
    
    Te_in  = args[0][0]
    tau = args[1]

    R   = args[3] 

    vshock = args[5]
    s      = args[6]

    R1 = 1
    R2 = R

    omega0 = sample_spherical(r = 1)
    # omega0 = np.array([1,0,0])
    
    
    segs = np.linspace(R1+1e-6,R2,1000)  
    segs_prob = 1/((segs)**4)
    r =   choice(segs ,1,p=segs_prob/np.sum(segs_prob))[0]

    r0 = sample_spherical(r = r)

    
    rMag = abs(np.sqrt(sum([i**2 for i in r0])))
    
    Te_r = Te_func(rMag,Te_in)

    V = -1* 1000 * vshock * r0
    
    lam_random = velocity2lambda(V[0],lam_0)

    nu0 = lam_2_nu(lam_random)
    


    w0 = 1
    scatters = 0
    
    Ne0 = Ne(tau,R = 1)
    
    Ne_density = lambda x: Ne0 *(1/(1+x))**s
    
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
            
                
            Te_r = Te_func(rMag,Te_in)
            
            # Find the photons distance to boundary
            # l0 = BoundaryDistance(r0,omega0,R2)
            
            # Mean free path for travelling photon
            mfp = MFP(nu0,Ne =  Ne_density(rMag) ,Te = Te_r)
                
            # Li = EscapeProbability(l0, mfp)
            Li = np.exp(-1 * Ne_density(rMag)*sigma_t)
            
            # weighting
            w1  = w0 - w0*Li
            
            #  Get the MFP
            lambda_i = -1 * mfp * np.log(1 - np.random.random()*(1 - Li))

            # Momentum of Electron
            vMag = get_electron_velocity(Te_r)
            p = vMag*me
            eta = p/(me*c)
        
            gamma = np.sqrt((eta**2) + 1)
            # print(gamma)
        
            v_div_c = eta / gamma
            mu0 = np.dot(v_e,omega0)
        
        
            x = 2 * gamma * ((h * nu0) / (me*c**2)) * (1 - mu0*v_div_c )
        
            while True:
                
                #  Get posible direction of the scattering and see if 
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
                    break
        
            nu1 = (1/h) * (x1) *  (1 / (2*gamma*(1 - (mu1*(v_div_c))))) * me * (c**2)
            
            
            if np.isnan(nu1):
               
                raise Exception()
                

            if np.random.random() < 0.375*sigma_hat(x)*(1 - mu0*v_div_c): break
                
            
        pathTaken.append((r0[0],r0[1],r0[2]))
        
        r0 = r0 + lambda_i*omega0
        
        
        rMag = abs(np.sqrt(sum([i**2 for i in r0])))
        
        if np.isnan(rMag) :
            
            #  Photon has been absorbed - don't care about uit 
            photonNu = [np.nan]
            photonWeights = [np.nan]
            scatters = np.nan
            pathTaken = [np.nan]
            
            break
        
        if w1<epsilion: 
            
            # ignore photons travelling in the negative x direction 
            if ( r0[0] < 0) or scatters < 1:
                photonNu = [np.nan]
                photonWeights = [np.nan]
                scatters = np.nan
                pathTaken = [np.nan]
                
            break

        photonWeights.append(w0*Li)
        photonNu.append(nu0)
        
        

        scatters+=1
   
        
        w0 = w1
        omega0 = omega1
        nu0 = nu1
        
       
        

        

    return scatters, photonNu,photonWeights,pathTaken
       
    
    

# =============================================================================
#  Run with multiprocessing handler
# =============================================================================

def mp_handler(Te=10000,
               tau=10,
               tau_CSM = None,
               R = 100,
               vwind = 200,
               vshock = 2000,
               nPhotons = 1e6,
               s = 2,
               returnPath = False,
               randomselection = False,
               photonHist = None):
    

    lams = []
    scatters = []
    
    start_time = time.time()
    
    nCPU = mp.cpu_count()
    nCPU = 8
    
    chunksize, extra = divmod(int(nPhotons) , 4 *  nCPU)
    print('Running model on %d cpus\n' % nCPU)
    if extra:
        chunksize += 1
        
    args = [[Te,tau,tau_CSM,R,vwind,vshock,s,returnPath,randomselection,photonHist]] * int(nPhotons)

    pathTaken = []
    photonWeights =[]

    
    with mp.Pool(processes =  nCPU) as pool:   

        for res in tqdm(pool.imap_unordered(func=scatterPhoton,
                                            iterable=args,
                                            # chunksize = chunksize
                                            ),
                        total = int(nPhotons)):
        
            scatters.append(res[0])
            lams.extend(res[1])
            photonWeights.extend(res[2])

            pathTaken.append(res[3])
                
                

        pool.close()
        pool.join()
        
        for p in pool._pool:
            if p.exitcode is None:
                p.terminate()
    
    print(f'\n\nCompleted model w1th: Te = {Te}K,tau = {tau}, R = {R}')
    print('\n\tTime taken: %.1fs' % (time.time()-start_time))
    return scatters,lams,photonWeights,pathTaken

    
# =============================================================================
# Get parameters ready
# =============================================================================



modelParams = list(product(tau_range,Te_range,R_range,vsh,vwind,s))
shuffle( modelParams )

args = modelParams[0]



print('Total number of model is %d' % len(modelParams))


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
        
        
        print('\n[%d / %d] Scattering %.1e photons' % (counter,len(modelParams),nPhotons))
        
                
        
    
        print(f'\n\tRunning model  through CDS: Te = {Te} K, tau = {tau}, R = {R}, Vsh = {vshock} km/s\n')

        
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
        
        # stop
        photonLams = (c/photonNu)/1e-10
        dl = (max(photonLams) - min(photonLams))/250
        
        photonLams_bins = np.arange(min(photonLams),max(photonLams),dl)
        
        photonVels = lambda2velocity(photonLams, lam_0/1e-10)
        photonVels_bins = np.arange(min(photonVels),max(photonVels),dl)
        
# =============================================================================
#         Save the line profile to file
# =============================================================================
        flx = []
        for i  in photonLams_bins:
            idx = abs( photonLams - i) < dl
            flx.append(np.sum(photonWeights[idx]))
            
        
        
        
        ax1.step(photonLams_bins- lam_0/1e-10,
                  flx,label = fr'$Tau={tau} || V = {vshock} || R = {R}$ ',lw = 0.5)
        
        ax1.step(-1*(photonLams_bins- lam_0/1e-10),flx,color = 'grey',ls = '--',lw = 0.25)



        print(f'\nRunning model Through CSM: Te = {Te} K, tau = {tau}, R = {R}, \n')
        

        
      
# =============================================================================
# 
# =============================================================================
        

        fpath = os.path.join(output_folder, fname)
              
        np.savetxt(fpath, np.vstack([photonLams_bins,flx/np.nanmax(flx)]).T,
                    delimiter=', ')
    
        print('\nSaved emission profile as: %s'% fpath)

    
# =============================================================================
    # plot velocity axis
# =============================================================================


        
        mn, mx = ax1.get_xlim()

           
        ax2.set_xlim(xmin = lambda2velocity(mn+lam_0/1e-10,lam0 = lam_0/1e-10),
                     xmax = lambda2velocity(mx+lam_0/1e-10,lam0 = lam_0/1e-10))
        

        
        ax2.axvline( -1 *vshock,lw = 0.5,ls = '--',color = 'red')
  
        
        trans = transforms.blended_transform_factory(ax2.transData, ax2.transAxes)
        
        ax2.text(-1 * vshock,0.01,'$V_{sh} = %.f~km/s$' % vshock,
                 rotation = 90,transform=trans,ha = 'right',
                 va = 'bottom',color = 'red',fontsize = 8)
        
        
        ax1.axvline(0 ,lw = 0.5,ls = ':',color = 'grey')
        ax1.legend(frameon = False,fontsize = 7)
        
                 
      
# =============================================================================
#      Histogram with number of scatters   
# =============================================================================

        scatters = np.array(scatters)[np.isfinite(scatters)]
        
        counts, bins = np.histogram(scatters )
        
        ax3.hist(scatters, histtype = 'step',lw = 0.5,color = 'black')
        
        ax3.axvline(tau**2,label ='$\\tau^2$',ls = '--',color = 'green',lw = 0.5)
        # ax3.axvline(tau,label ='$\\tau$',ls = ':',color = 'green',lw = 0.5)
        
        trans = transforms.blended_transform_factory(ax3.transData, ax3.transAxes)
        ax3.text(tau**2,0.01,'$\\tau_e^2$',rotation = 90,transform=trans,
                 ha = 'left',va = 'bottom',color = 'green')
        # ax3.text(tau,0.01,'$\\tau_e$',rotation = 90,transform=trans,ha = 'left',va = 'bottom',color = 'green')
        handles, labels = ax1.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        
# =============================================================================
#         Save the output profile
# =============================================================================

        
        ax1.set_xlabel('$\lambda~[A]$')
        ax2.set_xlabel('Velocity [$km~s^{-1}$]')
        ax3.set_xlabel('nScatters')
        
        
        # ax1.set_yscale('log')
        # ax1.set_ylim(1,ax1.get_ylim()[1])
        ax3.set_yscale('log')
        ax3.set_xscale('log')
        
        ax1.set_ylabel('$N$')
        ax3.set_ylabel('$Occurance$')
        
        fig.tight_layout()
        
        fpathFig = os.path.join(figures_folder, fr'tau_{tau}_T_e_{Te}_R_{R}_vs_{vshock}_vw_{vwind}_s_{s}.png')
        plt.savefig(fpathFig,dpi = 300)
        
                # plt.close('all')
                
                
        plt.show(block = False)
                # break
                            
                                 
                    