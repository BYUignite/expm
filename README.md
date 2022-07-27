# Exponential matrix in C++
* Based on [Julia's expm](https://github.com/JuliaLang/julia/blob/d386e40c17d43b79fc89d3e579fc04547241787c/base/linalg/dense.jl#L395-L422)
* Uses N. Higham "The scaling and squaring method of the matrix exponential revisited," SIAM Review [51(4):747-764, 2009](http://dx.doi.org/10.1137/090768539)

* Tested against Matlab and Python expm
* Requires Boost and LAPACK
* Compile Mac: g++ -std=c++11 -I/opt/homebrew/include -framework Accelerate expm.cc main.cc
    * the -I/opt/homebrew/include points to the boost header files
    * the -framework Accelerate links against LAPACK.

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


