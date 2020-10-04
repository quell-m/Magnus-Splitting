
include("rosen-zener.jl")

const testcase="Rosen-Zener"
const t0 = -5.0 # start time
const tend = 5.0 # end time
const dt = 0.1 # first step size
const reps=100 # number of iterations for benchmark

include("main.jl")
