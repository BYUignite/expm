# Exponential matrix in C++

Uses two similar methods for computing the matrix exponential.

1. [Julia's expm](https://github.com/JuliaLang/julia/blob/d386e40c17d43b79fc89d3e579fc04547241787c/base/linalg/dense.jl#L395-L422) (MIT license)
    * Nicholas J. Higham, "The squaring and scaling method for the matrix exponential revisited", SIAM Journal on Matrix Analysis and Applications, 26(4), 2005, 1179-1193.
        * Reprinted online SIAM Review Vol. [51(4):747-764](http://dx.doi.org/10.1137/090768539) (2009) 
2. Najfeld and Havel, Derivatives of the Matrix Exponential and Their Computation, Advances in Applied Mathematics, [16:321-375](https://doi.org/10.1006/aama.1995.1017) (1995).
    * This is referenced and discussed in Higham 2005, cited above.

Both of these use a scaling and squaring algorithm with Pade approximations.

* Tested against Matlab and Python ```expm```
* Requires Boost and LAPACK
* Compile Mac: ```g++ -std=c++11 -I/opt/homebrew/include -framework Accelerate expm.cc main.cc```
    * ```-I/opt/homebrew/include``` points to the boost header files
    * ```-framework Accelerate``` links against LAPACK.

## ODE Example

Consider vector $y(t)$ and constant matrix $A$

$$\frac{dy}{dt} = Ay,$$

$$y(0) = y_0.$$

Let $A$ be written in terms of its eigendecomposition

$$A = V\Lambda V^{-1}.$$

Then

$$\frac{dy}{dt} = V\Lambda V^{-1}y,$$

$$V^{-1}\frac{dy}{dt} = \Lambda V^{-1}y.$$

Let $\hat{y}=V^{-1}y$:

$$\frac{d\hat{y}}{dt} = \Lambda\hat{y}.$$

This has solution components

$$\hat{y}_i = \hat{y}_{i,0}e^{\Lambda_{i,i}t},$$

or as a vector,

$$\hat{y} = e^{\Lambda t}\hat{y}_0,$$

where $e^{\Lambda t}$ is diagonal with diagonal elements $e^{\Lambda_{i,i}t}$. Then

$$y = Ve^{\Lambda t}V^{-1}y_0,$$

$$y = e^{A t}y_0,$$

where 

$$e^{At} = Ve^{\Lambda t}V^{-1}.$$


## Note
* The boost matricies can by initialized from std::vectors, as, e.g., 
    ```A.data() = v;``` where ```A``` is a boost matrix and ```v``` is an ```std::vector```.
    * This assumes ```std::vector<foo>``` is used as the type of the storage array (third parameter of the boost matrix template) is ```std::vector<foo>```.
    * Otherwise, ```std::copy(v.begin(), v.end(), A.data().begin())``` works too, where ```copy``` is from ```#include<algorithm>```.
