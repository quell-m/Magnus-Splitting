using LinearAlgebra
using SparseArrays
using FFTW

include("expmv.jl")

function call_expm!(dt::Float64, B, u::Array{Complex{Float64}})
    expmv!(-1im * dt, B, u, tol=1e-12, m=30)
    nothing
end

include("time_stepper.jl")
const k = 50
const x = collect(1:2 * k)

const omega = 0.5
const T_0 = 1.0
const V_0 = 1.0

function F1(t::Float64)
    (V_0 * cos(omega * t) / cosh(t / T_0))
end
function F1d(t::Float64)
    (-V_0 * omega * sin(omega * t) / cosh(t / T_0) - V_0 * cos(omega * t) * sinh(t / T_0) / (cosh(t / T_0)^2 * T_0))
end
function F2(t::Float64)
    (-V_0 * sin(omega * t) / cosh(t / T_0))
end
function F2d(t::Float64)
    (-V_0 * omega * cos(omega * t) / cosh(t / T_0) + V_0 * sin(omega * t) * sinh(t / T_0) / (cosh(t / T_0)^2 * T_0))
end
const S1 = [0.0 1.0 ;1.0 0.0]
const S2 = [0.0  -1.0im; 1.0im 0.0]              
const H1 = kron(S1, sparse(1.0I, k, k))
const H2 = kron(S2, spdiagm(-1 => ones(k - 1), 1 => ones(k - 1)))

struct RosenMatrix <: TimeDependentMatrix
    c_H1::Complex{Float64}
    c_H2::Complex{Float64}
end

function evaluate_A(t::Float64)
    return RosenMatrix(F1(t),F2(t))
end

function evaluate_A_t(t::Float64)
    return RosenMatrix(F1d(t),F2d(t))
end

function evaluate_A(t::Array{Float64}, c::Array{Float64})
    c_H1 = 0.0
    c_H2 = 0.0
    for i in 1:length(t)
        c_H1 += c[i] * F1(t[i])
        c_H2 += c[i] * F2(t[i])
    end
    return RosenMatrix(c_H1,c_H2) # multiply by -1im
end

function evaluate_A_t(t::Array{Float64}, c::Array{Float64})
    c_H1 = 0.0
    c_H2 = 0.0
    for i in 1:length(t)
        c_H1 += c[i] * F1d(t[i])
        c_H2 += c[i] * F2d(t[i])
    end
    return RosenMatrix(c_H1,c_H2) # multiply by -1im
end

function setPsi0(t0::Float64)
    Psi(ones(2 * k), t0)
end

function LinearAlgebra.mul!(Y, H::RosenMatrix, B)
    Y[:] = H.c_H1*(H1*B) + H.c_H2*(H2*B)
    global mulc
    mulc+=1
    nothing
end

function Base. *(A::Number, B::RosenMatrix) 
    return RosenMatrix(A*B.c_H1,A.*B.c_H2)
end

function Base. +(A::RosenMatrix, B::RosenMatrix) 
    return RosenMatrix(A.c_H1+B.c_H1,A.c_H2+B.c_H2)
end

function LinearAlgebra.size(H::RosenMatrix, dim::Int)
     dim < 1 ? error("arraysize: dimension out of range") :
   return  (dim < 3 ? 2*k : 1)
end

function LinearAlgebra.ishermitian(H::RosenMatrix)
    return true
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