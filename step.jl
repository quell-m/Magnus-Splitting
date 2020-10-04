include("psi.jl")
include("schemes.jl")

using LinearAlgebra

mulc = 0 # counts the number of matrxi vector multiplications

abstract type TimeDependentMatrix end
using Base

function Base. *(A::TimeDependentMatrix, B) 
    Y = similar(B)
    mul!(Y, A, B)
    return Y
end

function step!(psi::Psi, dt::Float64, scheme::SplittingScheme, operator_sequence="AB")
    for k = 1:length(scheme.scheme)
        if scheme.scheme[k] != 0.0
            which_operator = operator_sequence[mod(k - 1, length(operator_sequence)) + 1]
            if which_operator == 'A'
                propagate_A!(psi, dt * scheme.scheme[k]) 
            elseif which_operator == 'B'
                propagate_B!(psi, dt * scheme.scheme[k]) 
            else
                error()
            end
        end
    end 
end

function step!(psi::Psi, dt::Float64, scheme::MagnusScheme)
    for j in 1:size(scheme.a)[1] # J
        B = evaluate_A(psi.t .+ dt .* scheme.c[j,:], scheme.a[j,:])
        call_expm!(dt, B, psi.u)
    end
    psi.t += dt
end

# psi solution, psi2 errorestimator
function step_adaptiv!(psi::Psi, psi2::Psi, dt, scheme::PalindromicScheme) 
    assign!(psi2, psi)
    step!(psi, dt, scheme, scheme.operator_sequence)
    step!(psi2, dt, scheme, reverse(scheme.operator_sequence))
    nothing
end

function compute_err(psi::Psi, psi2::Psi, dt, tol, scheme::PalindromicScheme) 
    err = 0.5 * distance(psi, psi2) / tol
end

function step_adaptiv!(psi::Psi, psi2::Psi, dt, scheme::DefectBasedScheme) 
    set_0!(psi2)
    which_operator = 'X'
    os = scheme.operator_sequence
    for k = 1:length(scheme.scheme)
        which_operator = os[mod(k - 1, length(os)) + 1]
        y = scheme.scheme[k]
        dt_k = dt * y
        if (k == length(scheme.scheme))
            y -= 1.0
        end
        if (which_operator == 'A')
            if (y != 0.0)
                add_apply_A!(psi, psi2, y)
            end
            if (dt_k != 0.0)
                propagate_A_derivative!(psi, psi2, dt_k)
            end
        elseif (which_operator == 'B')
            if (dt_k != 0.0)
                propagate_B_derivative!(psi, psi2, dt_k, scheme.order, k)
            end
            if (y != 0.0)
                add_apply_B!(psi, psi2, y)
            end
        else
            error()
        end 

    end
    if 'A' in os && which_operator != 'A' 
        add_apply_A!(psi, psi2, -1.0)
    end
    if 'B' in os && which_operator != 'B' 
        add_apply_B!(psi, psi2, -1.0)
    end
    nothing
end

function compute_err(psi::Psi, psi2::Psi, dt, tol, scheme::DefectBasedScheme) 
    err = dt * norm(psi2) / (scheme.order + 1) / tol
end

function step_adaptiv!(psi::Psi, psi2::Psi, dt, scheme::EmbeddedScheme) 
    kk = -1
    for k = 1:length(scheme.scheme)
        if (kk < 0 && k <= length(scheme.scheme2) && scheme.scheme[k] != scheme.scheme2[k])
            kk = k
            assign!(psi2, psi)
        end
        if (scheme.scheme[k] != 0.0)
            which_operator = scheme.operator_sequence[mod(k - 1, length(scheme.operator_sequence)) + 1]
            if which_operator == 'A'
                propagate_A!(psi, dt * scheme.scheme[k]) 
            elseif which_operator == 'B'
                propagate_B!(psi, dt * scheme.scheme[k]) 
            else
                error()
            end
        end
    end 
    if (kk < 0)
        kk = length(scheme.scheme) + 1
        assign!(psi2, psi)
    end   
    for k = kk:length(scheme.scheme2)
        if (scheme.scheme2[k] != 0.0)
            which_operator = scheme.operator_sequence[mod(k - 1, length(scheme.operator_sequence)) + 1]
            if which_operator == 'A'
                propagate_A!(psi2, dt * scheme.scheme2[k])
            elseif which_operator == 'B'
                propagate_B!(psi2, dt * scheme.scheme2[k])
            else
                error()
            end
        end 
    end     
    nothing
end
function compute_err(psi::Psi, psi2::Psi, dt, tol, scheme::EmbeddedScheme) 
    err = distance(psi, psi2) / tol
end
function step_adaptiv!(psi::Psi, psi2::Psi, dt, scheme::MilneScheme) 
    kk = -1
    for k = 1:length(scheme.scheme)
        if (kk < 0 && k <= length(scheme.scheme2) && scheme.scheme[k] != scheme.scheme2[k])
            kk = k
            assign!(psi2, psi)
        end
        if (scheme.scheme[k] != 0.0)
            which_operator = scheme.operator_sequence[mod(k - 1, length(scheme.operator_sequence)) + 1]
            if which_operator == 'A'
                propagate_A!(psi, dt * scheme.scheme[k]) 
            elseif which_operator == 'B'
                propagate_B!(psi, dt * scheme.scheme[k]) 
            else
                error()
            end
        end
    end 
    if (kk < 0)
        kk = length(scheme.scheme) + 1
        assign!(psi2, psi)
    end   
    for k = kk:length(scheme.scheme2)
        if (scheme.scheme2[k] != 0.0)
            which_operator = scheme.operator_sequence[mod(k - 1, length(scheme.operator_sequence)) + 1]
            if which_operator == 'A'
                propagate_A!(psi2, dt * scheme.scheme2[k])
            elseif which_operator == 'B'
                propagate_B!(psi2, dt * scheme.scheme2[k])
            else
                error()
            end
        end 
    end   
    nothing  
end
function compute_err(psi::Psi, psi2::Psi, dt, tol, scheme::MilneScheme) 
    err = scheme.k * distance(psi, psi2) / tol
end

function sum_ad(o, x, p::Integer, b::Array{Complex{Float64},1})
    xb = x * b
    ss = xb
    for j in 2:p
        s = xb
        ob = b
        for i in 1:j - 1
            ob = o * ob
            s = o * s + (binomial(j - 1, i) * (-1)^i) * x * ob
        end
        ss +=  1 / (factorial(j)) * s
    end   
    return ss
end

function step_adaptiv!(psi::Psi, psi2::Psi, dt, scheme::MagnusScheme) 
    B = evaluate_A(psi.t .+ dt .* scheme.c[1,:], scheme.a[1,:])# first element j
    B_t = evaluate_A_t(psi.t .+ dt .* scheme.c[1,:], scheme.c[1,:] .* scheme.a[1,:])
    call_expm!(dt, B, psi.u) # step
    psi2.u[:] = sum_ad(-1im * dt * B, -1im * B + -1im * dt * B_t, scheme.order, psi.u)
    for j in 2:size(scheme.a)[1]
        B = evaluate_A(psi.t .+ dt .* scheme.c[j,:], scheme.a[j,:])
        B_t = evaluate_A_t(psi.t .+ dt .* scheme.c[j,:], scheme.c[j,:] .* scheme.a[j,:])
        call_expm!(dt, B, psi2.u)
        call_expm!(dt, B, psi.u) # step
        psi2.u[:] += sum_ad(-1im * dt * B, -1im * B + -1im * dt * B_t, scheme.order, psi.u)
    end
    psi2.t = 0.0
    psi2.u[:] -= -1im * evaluate_A(psi.t + dt) * psi.u
    psi.t = psi.t + dt
    nothing
end

function compute_err(psi::Psi, psi2::Psi, dt, tol, scheme::MagnusScheme) 
    err = dt * norm(psi2) / (scheme.order + 1) / tol
end
