# ML-Flux
ML-Flux is a machine-learning software for determination of metabolic flux distributions from isotope labeling patterns. ML-Flux is available for direct use online on metabolicflux.org, where there are no specific system requiremts or need for installation of software. A detailed manual for online-use, as well as a demo dataset, is available on the website. The code for the core software used on metabolicflux.org is available in this repository.

# System requirements and installation
The following Python dependencies (and the most recently used versions) are required for local operation:
* keras (2.10.0)
* keras-preprocessing (1.1.2)
* numpy (1.26.4)
* pandas (2.2.2)
* scikit-learn (1.5.1)
* scipy (1.14.0)
* tensorflow (2.10.0)
There are no non-standard hardware dependencies.

There are no specific installation instructions. The code may be directly cloned from this repository

# Instructions for use and demo
Flux predictions may be conducted from runML-Flux.py in the Flux_Prediction folder. This file is dependent on other files within the entire repository and draws from models in the Trained_Models folder. Demo data is available on metabolicflux.org, or scripts for generating new test/demo data are provided in Data_Generation_and_Simulation



