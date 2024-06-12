

[arxiv_link]: https://arxiv.org/abs/2312.13280


[arxiv_link_SN2023fyq]: https://arxiv.org/abs/2401.15148




# eScatter.py - Electron Scattering in Python

*eScatter.py* is a Monte Carlo Electron Scattering Simulation designed to perform Monte Carlo simulations of electron scattering events.

This code was developed in order to better understand the emission lines from the interacting supernova, [SN 2021adxl][arxiv_link] and [SN 2023fyq ][arxiv_link_SN2023fyq], specifically the blue excess seen in the  H&alpha; 6563A  and He I 5876 emission line.

Scattering equations are adapted from [Pozdnyakovet al. 1983](https://ui.adsabs.harvard.edu/abs/1983ASPRv...2..189P/abstract), with similar codes used for interacting supernovae such as [SN 2010jl](https://arxiv.org/abs/1312.6617) and [SN 2013L](https://arxiv.org/abs/2003.09709).

> [!NOTE]
> The code is under active development through colloborations with Harvard University - Please direct any questions to  [me](mailto:sean.brennan@astro.su.se)



## Toy model

This code was build for to model the interacting of fat moving ejecta produced in a Supernova explosion colliding with slow moving material (see [Section 4 in our paper][arxiv_link]).

In short this code will  follow a photon, which was formed in a thin interface between the supernova ejecta and surrounding material, as it travels radially outwards through the dense material (described by an given optical depth, &tau;). As the photon travels outwards through this material, it will scatter of a electron with a certain energy (related to its Maxwellian velocity). The photon will the  undergo probabilistic iterations and may (or may not) scatter of the electron.



<p align="center">
  <img src="./eScatter_6563_models.png" alt="Image">
  <br>
  <em>Figure 1: An example of *eScatter.py*'s output for the H&alpha; emission line. Here an input profile similar to a thin emitting shell travelling at 3000 km/s is used. Photons are allowed to travel through a medium with a varying optical depth. Generally, with higher optical depth, the emergent profile is broader and the peak moves towards the blue. See the above link to see how we fit this to an observed  H&alpha; profile.</em>
</p>




The photon will continue to scatter outwards until it reaches an optically thin region, far away from where it was formed. We assume that photons emitted (or more accurately the photons which escape) on a single hemisphere would are observed (i.e. photons which escape and are travelling away from the observer are not seen). A histogram of the emergent photons is then take to be the emergent spectral profile.


> [!CAUTION]
> Understand the output of the code before applying models to data

Although scattering can produce very broad profiles when you assume high optical depths, you have to account for diffusion times. i.e. has a transient evolved for long enough for such photons to diffusion to optical thin regions. Typically Type IIn SNe will show such features post peak (after a few weeks post explosion), constaining the optical depth.

The diffusion time for a photon traveling through a medium with a certain optical depth, and assuming a random walk, is given by the following equation:

 $$ t_d = \frac{R^2 \tau}{c} $$

 where:
- \( t_d \) is the diffusion time,
- \( R \) is the characteristic length scale of the medium,
- \( \tau \) is the optical depth of the medium,
- \( c \) is the speed of light in a vacuum.








## Usage

This version of escatter.py is executed by changing values in the script and running the file. Below is the snippet of the parameters, many of these can be made to list to build a grid of models.


```python
# Rest wavelngth in units of Angstrom
lam_0 = 6563

# Denisty parameter either 0 or 2 [0 untested]
s = [2]

# Number of photons to send out
nPhotons = 1e5

# Optical depth (can be list)
tau_range = [3,5,10,15,25]

# Electron temperature from just above formation site to edge of boundary
Te_range = [(2e4,1e4)]

# Inner and outer radius of scattering region (this doesn't change the profile too much)
R_range = [(1e14,1.01e14)]

# wind velocity (untested)
vwind = [40] #km/s

# Shock velocity in km/s
vsh = [3000]

```

The above paramters produced the output used in Figure. 1. Once you have input your parameters, run the scipt using:

```bash
python escatter.py
```

Figures and output txt files will be written to `output/` and `figures/`

## Limitations

There are many limitations and caveats to using this code:

- The code takes no account of physical distances i.e. knowing where the shock from and how far the CSM extends is non-trivial.
- The code is completely spherically symmetric - such asymmetric features such as torus or jets, in reality may yield different results.
- Several functions are hardcoded, such as how and where the photons are injected, and how the CSM density and temperature evolves. This is a limitation of the code.
- The code does not include absorption (i.e. profile cannot produce P-Cygni-like profiles) - **under development**.

## Requirements
```
matplotlib==3.8.0
numpy==1.26.2
tqdm==4.66.1
```

## Warning
I make no claim as to the validity of the output of this code, use at you own discretion.
