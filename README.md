# Moment varieties from inverse Gaussian and gamma distributions
This repository contains auxiliary files for the manuscript _Moment varieties from inverse Gaussian and gamma distributions_ by Oskar Henriksson, Lisa Seccia, and Teresa Yu.

## Contents
The respository consists of the following directories:

- `ranks` with Maple files `rank_gamma.mpl` and `rank_inverse_gaussian.mpl` for the rank computations discussed in §5.1 of the paper. These ranks prove that the $k$-th secant variety of $\mathcal{M}_d$ has the expected dimension for $k,d\leqslant 100$. Explicit choices of parameters and the resulting Jacobians are given in the files `ranks_jacobians_gamma.txt` and `jacobians_inverse_gaussian.txt` for the cases with $d\leq 8$, since these are used as the base case for the paper [Moment varieties of the inverse Gaussian and gamma distributions are nondefective](https://arxiv.org/abs/2409.18421) by Oskar Henriksson, Kristian Ranestad, Lisa Seccia, and Teresa Yu.

- `identifiability_degrees` with Julia files `gamma_identifiability_degree.jl` and `inverse_gaussian_identifiability_degree.jl` for numerical estimation of the identifiability degrees, as well as a file `certify_orbits.jl` with a function for certifying orbits of a group action. Numerical certificates of the computed solutions are also included in the directory.

- `ed_degrees` with Julia files `edd_groebner.jl` and `edd_nag.jl` for symbolic and numeric computation, respectively, of the ED degrees given in §5.3 of the paper.

  
## Dependencies
The rank computations were done in Maple 2023. For the numerical degree compuations, HomotopyContinuation.jl v2.9.3 was used, whereas the symbolic degree computations were done in Oscar v0.13.0 (see `julia_environment/Manifest.toml` for details). 
