include("step.jl")
########################################################################
function equidistant_time_stepper(psi::Psi, t0::Float64, tend::Float64, dt::Float64, 
    scheme, out=false, outPhi=false)
    t = t0
    if (out)
        global counter
        counter += 1
    end
    vec_t = []
    vec_dt = []
    plotdata = abs.(psi.u)
    while t < tend
        curr_dt = min(dt, tend - t)
        step!(psi, curr_dt, scheme)
        if (out)
            push!(vec_t, t)
            push!(vec_dt, dt)
            if (outPhi)
                plotdata = [plotdata abs.(psi.u)]
            end
        end
        t += curr_dt
    end
    if (out)
        figure()
        plot(vec_t, vec_dt, label=scheme.name)
        savefig("time_$(counter).png")
        if (outPhi)
            figure()
            pcolormesh(push!(vec_t, tend), x, plotdata)
            savefig("spacetime_$(counter).png")
        end
    end
end

########################################################################

function adaptive_time_stepper(psi::Psi, t0::Float64, tend::Float64, dt::Float64,
tol::Float64, scheme, out=false, outPhi=false)
    psi2 = deepcopy(psi) 
    psi0 = deepcopy(psi) 
    t = t0
    facmin = 0.25
    facmax = 4.0
    fac = 0.9
    nsteps = 0
    nrejec = 0
    if (out)
        global counter
        counter += 1
    end
    vec_t = []
    vec_dt = []
    plotdata = abs.(psi.u)
    while (t < tend)
        nsteps += 1
        dt0 = dt
        err = 2.0
        while err >= 1.0
            dt = min(dt, tend - t)
            dt0 = dt
            step_adaptiv!(psi, psi2, dt, scheme)
            err = compute_err(psi, psi2, dt, tol, scheme)
            dt = dt * min(facmax, max(facmin, fac * (1.0 / err)^(1.0 / (scheme.order + 1))))
            if err >= 1.0
                assign!(psi, psi0)
                nrejec += 1
                # println("rejected: ", t, " ", dt, " ", err)          
            else
                if (out)
                    push!(vec_t, t)
                    push!(vec_dt, dt0)
                    if (outPhi)
                        plotdata = [plotdata abs.(psi.u)]
                    end
                end
                t += dt0
                assign!(psi0, psi)
            end
        end
    end
   # println("steps: ", nsteps," ",nrejec, ", Scheme: ", scheme.name)
    if (out)
        figure()
        plot(vec_t, vec_dt, label=scheme.name)
        savefig("time_$(counter).png")
        if (outPhi)
            figure()
            pcolormesh(push!(vec_t, tend), x, plotdata)
            savefig("spacetime_$(counter).png")
        end
    end
    [nsteps, nrejec]
end