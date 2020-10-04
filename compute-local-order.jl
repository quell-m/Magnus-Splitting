using Printf

# include("quantum-control.jl"); const ϵ = 2.0^(-8) # parameter
include("rosen-zener.jl")
# include("hubbard.jl")

function local_orders(psi0::Psi, t0::Float64, dt::Float64, 
                       scheme; ref_scheme=PP56A, operator_sequence="AB", rows=15)
    tab = zeros(Float64, rows, 5)
    nroot = 1
    psi = deepcopy(psi0)
    psi2 = deepcopy(psi0)
    psi_ref = deepcopy(psi0)
    dt1 = dt
    err_old1 = 0.0 
    err_old2 = 0.0
    println(" &           dt    &      err     &    p    &         err  &    p   \\\\ \\hline")
    @printf("%%---%s-----------------------------------------------------\n",scheme.name)
    for row = 1:rows
        if typeof(scheme) == PalindromicScheme
            step_adaptiv!(psi, psi2, dt1, scheme)
            assign!(psi2, 0.5 * (psi2 + psi))
        elseif typeof(scheme) == DefectBasedScheme
            step_adaptiv!(psi, psi2, dt1, scheme)
            assign!(psi2, psi - (dt1 / (scheme.order + 1)) * psi2)
        elseif typeof(scheme) == EmbeddedScheme
            step_adaptiv!(psi, psi2, dt1, scheme)
        elseif typeof(scheme) == MilneScheme
            step_adaptiv!(psi, psi2, dt1, scheme)
            assign!(psi2, psi - (scheme.k) * (psi2 - psi))
        elseif typeof(scheme) == MagnusScheme
            step_adaptiv!(psi, psi2, dt1, scheme)
            assign!(psi2, psi - (dt1 / (scheme.order + 1)) * psi2)
        else
            error("No suitable scheme type")
        end
        assign!(psi_ref, psi0)
        equidistant_time_stepper(psi_ref, t0, t0 + dt1, dt1 / 10, ref_scheme)
        err1 = distance(psi, psi_ref)
        err2 = distance(psi2, psi_ref)
        if (row == 1)
            @printf("%3i & %12.3e & %12.3e &         & %12.3e &        \\\\ \n", row, Float64(dt1), Float64(err1), Float64(err2))
            tab[row,1] = dt1
            tab[row,2] = err1
            tab[row,3] = 0 
            tab[row,4] = err2
            tab[row,5] = 0 
        else
            p1 = log(err_old1 / err1) / log((2.0)^(1.0 / nroot));
            p2 = log(err_old2 / err2) / log((2.0)^(1.0 / nroot));
            @printf("%3i & %12.3e & %12.3e & %7.2f & %12.3e &%7.2f \\\\ \n", row, Float64(dt1), Float64(err1), Float64(p1), Float64(err2), Float64(p2))
            tab[row,1] = dt1
            tab[row,2] = err1
            tab[row,3] = p1
            tab[row,4] = err2
            tab[row,5] = p2 
        end
        err_old1 = err1
        err_old2 = err2
        dt1 = 1 / (2.0)^(1.0 / nroot) * dt1
        assign!(psi, psi0)
    end
    println("\\hline\n")
    tab
end


dt = 1.0e-0
t0 = 0.01
schemes = [Magnus4, PP56A_D, PP34A_D, CF4oH, CF6n] 

psi0 = setPsi0(t0)
for scheme in schemes
    local_orders(psi0, t0, dt, scheme; ref_scheme=PP56A, operator_sequence="AB", rows=15);
end