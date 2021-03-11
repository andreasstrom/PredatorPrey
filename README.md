# PredatorPrey
A finite element solution of the predator-prey equation (the Lotke-Voltera system of ODEs) in MATLAB and [FEniCS](https://fenicsproject.org/olddocs/dolfin/latest/python/index.html)

## in MATLAB
The two versions of the solution implemented in matlab is for a siplified version of the equation. 
mkI is implemented in one dimension, omitting time dependent and non-linear terms, while mkII only assumes a constant predator population *v* = 1.

## in FEniCS
The solution implemented in Python  for the following equation:

<img src="https://user-images.githubusercontent.com/62762922/110867495-1fe2eb80-82c7-11eb-80c2-0e6d1576a227.png" data-canonical-src="https://user-images.githubusercontent.com/62762922/110867495-1fe2eb80-82c7-11eb-80c2-0e6d1576a227.png" width="500" height="100" />

where delta_1, delta_2, alpha, and beta are constants, f and g some given source terms, and u_0 and v_0 some given initial conditions.
