# Comparison Magnus und Splitting Integrators

Julia `1.4.2` has been used for the benchmarks.
You need the following additional packages

`import Pkg`

`Pkg.add("Statistics")`

`Pkg.add("FFTW)` 

`Pkg.add("JDL2")`


Run the benchmarks by executing the main files for the respective example

`julia main_hubbard.jl`

`julia main_quantum-control.jl`

`julia main_rosen-zener.jl`

The local error and error of the error estimator are computed by

`julia compute-local-order.jl`



