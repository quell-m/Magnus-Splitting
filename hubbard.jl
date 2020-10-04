using LinearAlgebra
import SparseArrays

include("expmv.jl")

function call_expm!(dt::Float64, B, u::Array{Complex{Float64}})
    expmv!(-1im * dt, B, u, tol=1e-12, m=30)
    nothing
end

include("time_stepper.jl")

#paramters for hubbard model
const x = collect(1:4900)
const ω = 3.5
const tp = 6.0
const σ = 2.0
const a = 0.2
const b = cos(ω * tp)

function F(t::Float64)
    exp((1im) * a * exp(-(t - tp)^2 / σ^2 / 0.2e1) *
        (cos(ω * (t - tp)) - cos(ω * tp)));
end
function Fd(t::Float64)
    ((-1im) * a * (t - tp) / (σ^2) * exp(-(t - tp)^2 / σ^2 / 0.2e1) *
        (cos(ω * (t - tp)) - cos(ω * tp)) +
        (-1im) * a * exp(-(t - tp)^2 / σ^2 / 0.2e1) *
        ω * sin(ω * (t - tp))) *
    exp((1im) * a * exp(-(t - tp)^2 / σ^2 / 0.2e1) * 
        (cos(ω * (t - tp)) - cos(ω * tp)))
end

using JLD2

@load "HHH.jld2" #matrizes H_symm, H_anti, and H_diag
@load "Groundstate.jld2" #initial data

println("Loaded static matrices")

struct HubbardMatrix <: TimeDependentMatrix
    c_symm::Complex{Float64}
    c_anti::Complex{Float64}
    c_diag::Complex{Float64}
end

function evaluate_A(t::Float64)
    return HubbardMatrix(real(F(t)),1im * imag(F(t)),1) # multiply by -1im
end

function evaluate_A_t(t::Float64)
    return HubbardMatrix(real(Fd(t)),1im * imag(Fd(t)),0) # multiply by -1im
end

function evaluate_A(t::Array{Float64}, c::Array{Float64})
    c_symm = 0.0
    c_anti = 0.0
    c_diag = 0.0
    for i in 1:length(t)
        c_symm += real(c[i] * F(t[i]))
        c_anti += imag(c[i] * F(t[i]))
        c_diag += c[i]
    end
    return HubbardMatrix(real(c_symm),1im * c_anti,c_diag) # multiply by -1im
end

function evaluate_A_t(t::Array{Float64}, c::Array{Float64})
    c_symm = 0.0
    c_anti = 0.0
    for i in 1:length(t)
        c_symm += real(c[i] * Fd(t[i]))
        c_anti += imag(c[i] * Fd(t[i]))
    end
    return HubbardMatrix(c_symm,1im * c_anti,0) # multiply by -1im
end

function setPsi0(t0::Float64=0.0)
    Psi(g, t0)
end

function LinearAlgebra.mul!(Y, H::HubbardMatrix, B)
    Y[:] = H.c_symm*(H_symm*B) + H.c_anti*(H_anti*B) + H.c_diag*(H_diag.*B)
    global mulc
    mulc+=1
    nothing
end

function Base. *(A::Number, B::HubbardMatrix) 
    return HubbardMatrix(A*B.c_symm,A*B.c_anti,A*B.c_diag)
end

function Base. +(A::HubbardMatrix, B::HubbardMatrix) 
    return HubbardMatrix(A.c_symm+B.c_symm,A.c_anti+B.c_anti,A.c_diag+B.c_diag)
end

function LinearAlgebra.size(H::HubbardMatrix, dim::Int)
     dim < 1 ? error("arraysize: dimension out of range") :
   return  (dim < 3 ? length(x) : 1)
end

function LinearAlgebra.ishermitian(H::HubbardMatrix)
    return imag(H.c_diag)==0&&imag(H.c_symm)==0&&real(H.c_anti)==0
end 

function propagate_A!(psi::Psi, dt::Float64)
    psi.t += dt
    nothing
end
function add_apply_A!(psi::Psi, psi2::Psi, dt::Float64)
    psi2.t += dt
    nothing
end
function propagate_A_derivative!(psi::Psi, psi2::Psi, dt::Float64)
    psi.t += dt
    nothing
end

function propagate_B!(psi::Psi, dt::Float64)
    call_expm!(dt, evaluate_A(psi.t), psi.u)
    nothing
end

function add_apply_B!(psi::Psi, psi2::Psi, dt::Float64)
    psi2.u[:] += -1im * dt * evaluate_A(psi.t) * psi.u
    nothing
end

function propagate_B_derivative!(psi::Psi, psi2::Psi, dt::Float64, order::Int64, k::Int)
    if (k > 2)
        call_expm!(dt, evaluate_A(psi.t), psi2.u)
    end
    call_expm!(dt, evaluate_A(psi.t), psi.u)
    psi2.u[:] += (sum_ad(-1im * dt * evaluate_A(psi.t), -1im * dt * evaluate_A_t(psi.t), order, psi.u) .* psi2.t)
    nothing
end