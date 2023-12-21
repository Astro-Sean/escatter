
[arxiv_link]: https://arxiv.org/abs/2312.13280




# Electron Scattering in Python

*escatter.py* is a Monte Carlo Electron Scattering Simulation designed to perform Monte Carlo simulations of electron scattering events.

This code was developed in order to better understand the emission lines from the interacting supernova SN 2021adxl [SN 2021adxl][arxiv_link], specifically the blue excess seen in the  H&alpha; 6563A emission line.



## Toy model

This code was build for one specific science case (see [Section 4 in our paper][arxiv_link]). In short this code will  follow a photon, which was formed in a thin interface between the supernova ejecta and surrounding material, as it travels radially outwards through the dense material (described by an given opacity, &tau;).

As the photon travels outwards through this material, it will scatter of a electron with a ceratin energy (related to its Maxwellian velocity). The photon will the  undergo a probablitic interation and may (or may not) scatter of the electron.

 The photon will continue to scatter outwards until it reaches an optically thin region, far away from where it was formed. We assume that photons emitted (or more accurately the photons which escape) on a single hemisphere would are observed (i.e. photons which escpae and are travelling away from the observer are not seen). A histogram of the emergent photons is then take to be the emergent spectral profile.

This code is adapted from [Pozdnyakovet al. 1983](https://ui.adsabs.harvard.edu/abs/1983ASPRv...2..189P/abstract), with similar codes being developed for interacting supernovae such as [SN 2010jl](https://arxiv.org/abs/1312.6617) and [SN 2013L](https://arxiv.org/abs/2003.09709).


## Usage

This version of escatter.py is executed by changing values in the script and running the file. Below is the snippet of the paramters, many of these can be made to list to build a grid of models.


```python
# rest wavelngth in units of Angstrom
lam_0 = 6563

# denisty parameter either 0 or 2 [0 untested]
s = [2]

# Number of photons to send out
nPhotons = 1e4

# Optical depth (can be list)
tau_range = [5]

# Electron temperature from just above formation site to edge of boundary, can be list
Te_range = [(30e3,10e3)]

#unitless radius of scattering region, cal be list
R_range = [10]

# wind velocity (untested), can be list
vwind = [0] #km/s

# Shock velocity in km/s, can be list
vsh = [5500]

```

Once you have input your parameters, run the scipt using:

```bash
python escatter.py
```

Figures and output txt files will be written to `output/` and `figures/`

## Limitations

There are many limitations and cavat to using this code:

- The code takes no account of physical distances i.e. knowing where the shock from and how far the CSM extends is non-trivial.
- The code is completely sphereically syemmetric - such assymetric featurs such as torus or jets, in reality may yield different results.
- Several functions are hardcoded, such as how and where the photons are injected, and how the CSM denisty and temperature evolves. This is a limitation of the code.

**If the community requires a more robust tool - reach out to [me](mailto:sean.brennan@astro.su.se)**
## Requirements
```
matplotlib==3.8.0
numpy==1.26.2
tqdm==4.66.1
```

## Warning
I make no claim as to the validity of the output of this code, use at you own descretion.
