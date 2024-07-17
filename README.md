# OperatorLearningforHPDEs

This repository holds the code used in the paper "Operator learning for hyperbolic partial differential equations," available on the [arXiv](https://arxiv.org/abs/2312.17489).

The figures used in the paper were generated with the default parameters in `main.m`. The code was run on an 64-core AMD Opteron with 512GB of RAM, using MATLAB version 9.13.0.2105380 (R2022b), Update 2.

## Dependencies

Our code requires [chebfun](https://www.chebfun.org/), which is available for download [here](https://www.chebfun.org/download/). Plotting some of the figures requires [inpaint_nans](https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans), an interpolation function written by John D'Errico.

## Usage 

Simply place `chebfun` and `inpaint\_nans.m` into the folder with the rest of the code. Add them to your path and run `main.m`. This runs Algorithm 2 of the aforementioned paper, saves all the relevant variables to a `.mat` file, and saves a figure containing three plots: 1) the approximate Green's function overlaid with partition blocks 2) the actual Green's function 3) the L2 error. *Note: The code will attempt to run in parallel; to avoid this, replace the function CONSTRUCTPAR with CONSTRUCT in `main.m`. Arguments will need to be modified.

Several parameters are available for modification, all of which are described in `main.m`. In theory, the Green's function `G` can be replaced by any anonymous function with the same arguments and output.

## Disclaimer

Our code is a proof-of-concept for a theoretical algorithm. It is not intended for real-world use for learning solution operators of hyperbolic PDEs.

Comments and suggestions welcome.
