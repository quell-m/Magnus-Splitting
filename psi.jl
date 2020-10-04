using LinearAlgebra

mutable struct Psi
    u::Array{Complex{Float64},1}
    t::Float64
end
function distance(p1::Psi, p2::Psi)
    sqrt(LinearAlgebra.norm(p1.u - p2.u)^2 + abs(p1.t - p2.t)^2)
end

function LinearAlgebra.norm(p1::Psi)
    sqrt(LinearAlgebra.norm(p1.u)^2 + abs(p1.t)^2)
end

function assign!(p1::Psi, p2::Psi)
    p1.u[:] = p2.u[:]
    p1.t = p2.t
    nothing
end

function set_0!(p1::Psi)
    p1.u .= 0.0+0.0im
    p1.t = 0.0
    nothing
end

import Base:+
function +(psi1::Psi,psi2::Psi)
    Psi(psi1.u+psi2.u,psi1.t+psi2.t)
end

import Base:*
function *(a::Number,psi1::Psi)
    Psi(a*psi1.u,a*psi1.t)
end

import Base:*
function *(psi1::Psi,a::Number)
    Psi(a*psi1.u,a*psi1.t)
end

import Base:-
function -(psi1::Psi,psi2::Psi)
    psi1+(-1)*psi2
end
