include("quantum-control.jl")

const testcase = "Quantum control"
const t0 = 0.0 # start time
const tend = 0.75 # end time
const dt = 0.01 # first step size
const reps = 100 # number of iterations for benchmark
#const ϵ = 2.0^(-6)
const ϵ = 2.0^(-8)
#const ϵ = 2.0^(-10)
#const ϵ = 2.0^(-12)

include("main.jl")