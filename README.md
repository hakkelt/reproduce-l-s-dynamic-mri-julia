# Efficient Dynamic Parallel MRI Reconstruction for the Low-Rank Plus Sparse Model

Julia implementation for reproducing the results in the paper: 
"Efficient Dynamic Parallel MRI Reconstruction for the Low-Rank Plus Sparse Model,"
IEEE Trans. on Computational Imaging, 5(1):17-26, 2019, by Claire Lin and Jeffrey A. Fessler,
EECS Department, University of Michigan

# Motivation

 - Make it possible to compare algorithms described in the paper to other algorithms implemented in Julia, and
 - Compare the following three Julia packages perfomance-wise: `LinearMaps.jl`, `AbstractOperators.jl`, `FunctionOperators.jl`
 
# Methods

 - Datasets:
   - **PINCAT**: L. Wissmann, C. Santelli, W. P. Segars, and S. Kozerke,
   “MRXCAT: Realistic numerical phantoms for cardiovascular magnetic resonance,”
   Journal of Cardiovascular Magnetic Resonance, vol. 16, no. 1, p. 63, Aug. 2014, doi: 10.1186/s12968-014-0063-3. 
   - **Multicoil cardiac cine MRI**, **Multicoil cardiac perfusion MRI**, **Multicoil abdominal dce MRI**:
   R. Otazo, E. Candès, and D. K. Sodickson, “Low-rank plus sparse matrix decomposition for accelerated dynamic MRI
   with separation of background and dynamic components,” Magnetic Resonance in Medicine, vol. 73, no. 3, pp. 1125–1136,
   2015, doi: 10.1002/mrm.25240. 
 - Algorithms:
   - **AL-CG**: Augmented Lagrangian with Conjugate Gradient step
   - **AL-2**: Optimized Augmented Lagrangian developed by Claire Lin and Jeffrey A. Fessler
   - **ISTA**: Proximal Gradient L+S (Iterative Shrinkage/Soft Thresholding Algorithm)
   - **FISTA**: Proximal Gradient L+S (Fast Iterative Shrinkage-Thresholding Algorithm) by Yurii Nesterov
   - **POGM**: Proximal Gradient L+S (Proximal Optimized Gradient Method), D. Kim and J. A. Fessler, “Optimized first-order methods for smooth convex minimization,” Math. Program., vol. 159, no. 1, pp. 81–107, Sep. 2016, doi: 10.1007/s10107-015-0949-3. 
 - Benchmarking:
   - In case of **Matlab code**, I measured only once using the `clock()` function (for some unknown reason,
   `tic`-`toc`, and `timeit()` gave wrong results)
   - In case of **Julia code**, I ran the algorithms 1+3 times (the first is the "warm-up"), and I measured the last three
   execution with BenchmarkTools.jl. At the results section below, I used the median value for running time, and I also added
   the aggregated memory usage of the algorithm measured also by BenchmarkTool.jl.

# Results

 - [MatlabToJulia.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/MatlabToJulia.ipynb) compares
 this Julia implementation to Lin and Fessler's reference implementation showing that they give technically identical results.
 
 - [MatlabOnly.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/MatlabOnly.ipynb)
 measures the speed of Lin and Fessler's reference implementation.
 
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 16.8 s | 49.0 s | 19.9 s | - |
  | AL-2 | 17.8 s | 55.3 s | 27.5 s | - |
  | ISTA | 1.7 s | 5.5 s | 2.3 s | 141.8 s |
  | FISTA | 2.4 s | 7.6 s | 3.1 s | 256.7 s |
  | POGM | 1.8 s | 5.8 s | 2.3 s | 140.0 s |

 - [LinearMaps_Naive.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/LinearMaps_Naive.ipynb)
 shows a "naive" implementation using LinearMaps. By "naive", I mean that no optimizations are applied.
 
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 30.6 s, 2.90 GiB | 79.4 s, 6.14 GiB | 33.3 s, 2.43 GiB | - |
  | AL-2 | 11.6 s, 12.20 GiB | 32.0 s, 33.17 GiB | 13.1 s, 13.82 GiB | - |
  | ISTA | 6.8 s, 3.40 GiB | 18.1 s, 8.67 GiB | 7.6 s, 3.62 GiB | 72.7 s, 6.71 GiB |
  | FISTA | 6.9 s, 3.40 GiB | 17.6 s, 8.67 GiB | 7.5 s, 3.62 GiB | 74.9 s, 6.71 GiB |
  | POGM | 6.9 s, 3.45 GiB | 17.8 s, 8.77 GiB | 7.5 s, 3.65 GiB | 73.7 s, 6.96 GiB |

 - [LinearMaps_Optimized.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/LinearMaps_Optimized.ipynb)
 diplays an optimized variant using in-place operations wherever it is possible reducing the number of allocations.
  
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 28.4 s, 11.26 GiB | 80.8 s, 31.08 GiB | 33.8 s, 12.95 GiB | - |
  | AL-2 | 9.8 s, 6.22 GiB | 28.4 s, 17.56 GiB | 11.0 s, 7.32 GiB | - |
  | ISTA | 6.0 s, 627.71 MiB | 19.6 s, 1.36 GiB | 6.8 s, 582.06 MiB | 72.7 s, 2.18 GiB |
  | FISTA | 6.2 s, 627.71 MiB | 16.5 s, 1.36 GiB | 6.8 s, 582.06 MiB | 73.2 s, 2.18 GiB |
  | POGM | 6.3 s, 677.71 MiB | 16.8 s, 1.45 GiB | 6.8 s, 622.06 MiB | 72.6 s, 2.42 GiB |
  
  - [AbstractOperators.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/AbstractOperators.ipynb)
 is a variant with AbstractOperators.jl that uses the same optimization tricks as LinearMaps_Optimized.ipynb. On the other
 hand, the code in this notebook is much more readable as the abundance of `vec` and `reshape` calls are reduced to minimal
 (LinearMaps.jl accepts only vectors as input, so I was forced to vectorize and de-vectorize matrices all the time I applied
 a LinearMap on them).
 
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 28.7 s, 678.06 MiB | 76.1 s, 1.45 GiB | 32.1 s, 622.41 MiB | - |
  | AL-2 | 9.0 s, 1.14 GiB | 23.3 s, 2.93 GiB | 10.2 s, 1.22 GiB | - |
  | ISTA | 6.4 s, 627.68 MiB | 16.4 s, 1.36 GiB | 7.3 s, 582.03 MiB | 72.8 s, 2.18 GiB |
  | FISTA | 6.5 s, 627.68 MiB | 16.3 s, 1.36 GiB | 7.3 s, 582.03 MiB | 72.0 s, 2.18 GiB |
  | POGM | 6.6 s, 677.68 MiB | 16.5 s, 1.45 GiB | 7.0 s, 622.03 MiB | 73.9 s, 2.42 GiB |
 
  - [FunctionOperators_Naive.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/FunctionOperators_Naive.ipynb)
  is a variant using FunctionOperators.jl. As the name suggests, it is similar to LinearMaps_Naive.ipynb, but it is "less
  naive" in the sense that the operators are optimized, and the algorithms are the unoptimized ("naive") parts of the notebook.
  This variant is the most readable meaning that it is the closest to the mathematical notation.
   
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 30.6 s, 2.90 GiB | 79.4 s, 6.14 GiB | 33.3 s, 2.43 GiB | - |
  | AL-2 | 11.6 s, 12.20 GiB | 32.0 s, 33.17 GiB | 13.1 s, 13.82 GiB | - |
  | ISTA | 6.8 s, 3.40 GiB | 18.1 s, 8.67 GiB | 7.6 s, 3.62 GiB | 72.7 s, 6.71 GiB |
  | FISTA | 6.9 s, 3.40 GiB | 17.6 s, 8.67 GiB | 7.5 s, 3.62 GiB | 74.9 s, 6.71 GiB |
  | POGM | 6.9 s, 3.45 GiB | 17.8 s, 8.77 GiB | 7.5 s, 3.65 GiB | 73.741 s, 6.96 GiB |
  
  - [FunctionOperators_Optimized.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/FunctionOperators_Optimized.ipynb)
  is similar AbstractOperators.ipynb as includes all possible optimization using in-place operations.
  
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 27.0 s, 640.69 MiB | 76.2 s, 1.38 GiB | 33.0 s, 592.54 MiB | - |
  | AL-2 | 8.6 s, 1.04 GiB | 22.7 s, 2.65 GiB | 10.0 s, 1.11 GiB | - |
  | ISTA | 6.1 s, 627.79 MiB | 16.6 s, 1.36 GiB | 7.6 s, 582.14 MiB | 75.1 s, 2.18 GiB |
  | FISTA | 6.1 s, 627.79 MiB | 16.5 s, 1.36 GiB | 7.9 s, 582.14 MiB | 73.8 s, 2.18 GiB |
  | POGM | 6.3 s, 677.79 MiB | 16.747 s, 1.45 GiB | 7.079 s, 622.14 MiB | 73.6 s, 2.42 GiB |
  
    - [FunctionOperators_Pretty.ipynb](https://github.com/hakkelt/reproduce-l-s-dynamic-mri-julia/blob/master/FunctionOperators_Pretty.ipynb)
  this variant is a mixture of FunctionOperators_Naive.ipynb and FunctionOperators_Optimized.ipynb: it combines the
  readability of FunctionOperators_Naive.ipynb and the memory-effectiveness of FunctionOperators_Optimized.ipynb by using
  the `@recycle` macro from FunctionOperators.jl.
  
  | algorithm \ data set | PINCAT | Multicoil cardiac cine MRI | Multicoil cardiac perfusion MRI | Multicoil abdominal dce MRI |
  | --- | --- | --- | --- | --- |
  | AL-CG | 28.0 s, 741.02 MiB | 75.1 s, 1.57 GiB | 30.7 s, 662.87 MiB | - |
  | AL-2 | 8.4 s, 1.26 GiB | 21.8 s, 3.26 GiB | 9.3 s, 1.36 GiB | - |
  | ISTA | 6.3 s, 1.05 GiB | 16.4 s, 2.49 GiB | 6.9 s, 1.04 GiB| 73.5 s, 3.03 GiB |
  | FISTA | 6.3 s, 1.05 GiB | 16.3 s, 2.49 GiB | 6.9 s, 1.04 GiB | 72.2 s, 3.03 GiB |
  | POGM | 6.4 s, 1.10 GiB | 16.5 s, 2.58 GiB | 6.9 s, 1.08 GiB | 74.5 s, 3.28 GiB |
  
  # Conclusion
  
   - The memory consumption of the packages are the following: `LinearMaps.jl` > `AbstractOperators.jl` > `FunctionOperators.jl`.
     - `LinearMaps.jl` is particularly inefficient memory-wise when non-square operators (i.e. operators with different input and
   output size) are composed as the composit operator re-allocates the buffer each time it is applied to a matrix.
     - `FunctionOperators.jl` is only slightly more memory-efficient as `AbstractOperators.jl` due to the fact that it reuses
     buffers when it is possible.
   - While I expected that the Julia code would be significatly faster than the Matlab, I succeded implementing only the
   AL-2 algorithm to be much faster than the Matlab implementation. Also, in case of Multicoil abdominal dce MRI data set,
   the Julia implementation is much faster because the Julia implementation of NFFT is much faster.
     - Later I might spend some time further optimizing my code to outperform the Matlab version.
     - It would be also interesting to compare the memory usage of Matlab functions, but that problem appears to be quite
     challenging (https://www.mathworks.com/matlabcentral/answers/97560-how-can-i-monitor-how-much-memory-matlab-is-using).
  
