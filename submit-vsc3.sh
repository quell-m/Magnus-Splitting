#!/bin/bash
#
#SBATCH -J julia
#SBATCH -N 1 # request an exclusive node
#SBATCH --mail-type=ALL    # first have to state the type of event to occur 
#SBATCH --mail-user=<quell@iue.tuwien.ac.at>   # and then your email address
#SBATCH --time=24:00:00

../julia-1.4.2/bin/julia main_rosen-zener.jl > rosen-zener.out
../julia-1.4.2/bin/julia main_hubbard.jl > hubbard.out
../julia-1.4.2/bin/julia main_quantum-control.jl > quantum-control.out
