import PyPlot
import SparseArrays
using FFTW
using LinearAlgebra

include("expmv.jl")

function call_expm!(dt::Float64, B, u)
    expmv!(-1.0im * dt, B, u, tol=1e-12)
    nothing
end

include("time_stepper.jl")
const nx = 2^11
const xmin = -1.0
const xmax = +1.0
roh(x) =   exp.(-1 ./ (1 .- x.^2)) .* (abs.(x) .< 1)
roh_t(x) =  .-(2.0 .* exp.(1 ./ (-1 .+ x.^2)) .* x) ./ ((-1 .+ x).^2 .* (1 .+ x).^2) .* (abs.(x) .< 1)

u0(x) = (1e-3 .* pi).^(-1 ./ 4) .* exp.(100 .* 1im .* (x .+ 0.3) .- 500 .* (x .+ 0.3).^2)
V0(x) = roh.(4 .* x) .* sin.(20 .* pi .* x)

VE(x,t) = V0.(x) .+ roh.(3 .* t .- 1) .* roh.(sin.(2 .* pi .* (x .- t))) 

VE_t(x,t) = 3 .* roh_t.(3 .* t .- 1) .* roh.(sin.(2 .* pi .* (x .- t))) .+ roh.(3 .* t .- 1) .* (-2 * pi) .* cos.(2 .* pi .* (x .- t)) .* roh_t.(sin.(2 .* pi .* (x .- t))) 

const x = collect(0:nx - 1) .* ((xmax - xmin) / nx) .+ xmin

kx = [0:nx / 2 - 1; 0.0; -nx / 2 + 1:-1]
kx = 1.0im * 2 * (pi / (xmax - xmin)) * kx
const kx2 = -kx .* kx # spectral coefficients

struct SemiclassicalMatrix <: TimeDependentMatrix
    off::Complex{Float64}
    diag::Array{Complex{Float64}}
end

function evaluate_A(t::Float64)
    return SemiclassicalMatrix(1,1.0 ./ ϵ .* VE.(x, t))# multiply by -1im
end
function evaluate_A_t(t::Float64)
    return SemiclassicalMatrix(0,1.0 ./ ϵ .* VE_t.(x, t))# multiply by -1im
end

function evaluate_A(t::Array{Float64}, c::Array{Float64})
    c_LL = 0.0
    c_diag = zeros(length(x))
    for i in 1:length(t)
        c_LL += c[i]
        c_diag[:] += c[i] * (1.0 ./ ϵ .* VE.(x, t[i]))
    end
    return SemiclassicalMatrix(c_LL,c_diag)
end

function evaluate_A_t(t::Array{Float64}, c::Array{Float64})
    c_diag = zeros(length(x))
    for i in 1:length(t)
        c_diag[:] += c[i] * (1.0 ./ ϵ .* VE_t.(x, t[i]))
    end
    return SemiclassicalMatrix(0,c_diag)
end

function setPsi0(t0::Float64=0.0)
    Psi(u0(x), t0)
end

psi = setPsi0()
#create plans for the FFT
const p_fft = plan_fft!(psi.u)
const p_ifft = plan_ifft!(psi.u)

function LinearAlgebra.mul!(Y, H::SemiclassicalMatrix, B)
    Y[:]=B[:]
    p_fft * Y
    Y .*= H.off * ϵ * kx2
    p_ifft * Y
    Y[:] += H.diag.*B
    global mulc
    mulc+=1
    nothing
end

function Base. *(A::Number, B::SemiclassicalMatrix) 
    return SemiclassicalMatrix(A*B.off,A.*B.diag)
end

function Base. +(A::SemiclassicalMatrix, B::SemiclassicalMatrix) 
    return SemiclassicalMatrix(A.off+B.off,A.diag.+B.diag)
end

function LinearAlgebra.size(H::SemiclassicalMatrix, dim::Int)
     dim < 1 ? error("arraysize: dimension out of range") :
   return  (dim < 3 ? nx : 1)
end

function LinearAlgebra.ishermitian(H::SemiclassicalMatrix)
    return imag(H.diag) == zeros(length(H.diag))
end 

function propagate_A!(psi::Psi, dt::Float64)
    p_fft * psi.u
    psi.u .*= exp.(dt .* -1.0im .* ϵ .* kx2)
    p_ifft * psi.u
    psi.t += dt
    global mulc
    mulc+=1
    nothing
end
function add_apply_A!(psi::Psi, psi2::Psi, dt::Float64)
    tmp = deepcopy(psi.u)
    p_fft * tmp
    tmp[:] .*= dt .* -1.0im .* ϵ .* kx2
    p_ifft * tmp
    psi2.u += tmp
    psi2.t += dt
    global mulc
    mulc+=1
    nothing
end
function propagate_A_derivative!(psi::Psi, psi2::Psi, dt::Float64)
    p_fft * psi.u
    p_fft * psi2.u
    psi2.u[:] .*= exp.(dt .* -1.0im .* ϵ .* kx2)
    psi.u[:] .*= exp.(dt .* -1.0im .* ϵ .* kx2)
    p_ifft * psi.u
    p_ifft * psi2.u
    global mulc
    mulc+=2
    psi.t += dt
    nothing
end

function propagate_B!(psi::Psi, dt::Float64)
    psi.u[:] .*= exp.(dt .* -1.0im ./ ϵ .* VE.(x, psi.t))
    nothing
end
function add_apply_B!(psi::Psi, psi2::Psi, dt::Float64)
    psi2.u[:] .+= ((dt .* -1.0im ./ ϵ .* VE.(x, psi.t)) .* psi.u)
    nothing
end
function propagate_B_derivative!(psi::Psi, psi2::Psi, dt::Float64, order::Int64, k::Int)
    psi.u[:] .*= exp.(dt .* -1.0im ./ ϵ .* VE.(x, psi.t))
    psi2.u[:] .= psi2.u .* exp.(dt .* -1.0im ./ ϵ .* VE.(x, psi.t)) .+ 
                 psi2.t .* (dt .* -1.0im ./ ϵ .* VE_t.(x, psi.t)) .* psi.u
    nothing
end