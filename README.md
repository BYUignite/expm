# Exponential matrix in C++
* Based on [Julia's expm](https://github.com/JuliaLang/julia/blob/d386e40c17d43b79fc89d3e579fc04547241787c/base/linalg/dense.jl#L395-L422)
* Uses N. Higham "The scaling and squaring method of the matrix exponential revisited," SIAM Review [51(4):747-764, 2009](http://dx.doi.org/10.1137/090768539)

* Tested against Matlab and Python expm
* Requires Boost and LAPACK
* Compile Mac: g++ -std=c++11 -I/opt/homebrew/include -framework Accelerate expm.cc main.cc

