# Electron Scattering Python

[arxiv_link]: https://arxiv.org/list/astro-ph/new

escatter is a Monte Carlo Electron Scattering Simulation designed to perform Monte Carlo simulations of electron scattering events.

This code was developed in order to better understand the emission lines from the interacting supernova SN2021adxl [SN 2021adxl][arxiv_link], specifically the blue excess seen in the  H\alpha 6563A emission line

## Limitations of this model

This code was build for one specific science case (see [Section 4 in our paper][arxiv_link] ). In short this code will  follow a photon, which was formed from a thin interface between the supernova ejecta and csm, radially outwards through dense matterial (described by an given opacity &tau: )
