include("hubbard.jl")

const testcase="Hubbard"
const t0 = 0.0 # start time
const tend = 30.0 # end time
const dt = 0.1 # first step size
const reps=100 # number of iterations for benchmark

include("main.jl")
