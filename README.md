# Torrey-Bloch_CubeSolver

Solving the Torrey–Bloch equation in a cubic domain.

This repository contains the code accompanying the research article:

> "Systematic Errors from Inhomogeneous Alkali Polarization in NMR Gyroscopes and Comagnetometers"

The main package is located in the `isotopeShift_Simu` folder, while the `figs` folder contains the code for calculating and plotting the figures presented in the article.

## Running the Codes

To run the codes, add the following folders to your MATLAB path:

- `isotopeShift_Simu`
- `isotopeShift_Simu/chebfun-master`
- `isotopeShift_Simu/specialFuns`

For demonstrations, refer to the following example scripts:

- `isotopeShift_Simu/+Table0/CubeSolver2_demo.m`
- `isotopeShift_Simu/+Table0/CubeSolver2_demo2_sweepCellTemperature_v2.m`

### Recommended MATLAB Version

The recommended MATLAB version is R2024a. Other versions newer than R2019b should work as well, though we have not performed comprehensive testing.