using Printf
using Statistics

println("Setup finished")

function bm(t0::Float64, tend::Float64, dt::Float64, reps::Int, testcase::String)
    schemes = [ EMB43AKp,EMB43AKp_D3,PP56A,PP56A_D,EMB54AKii,CF4oH,CF6n] 
    for tol in [1e-5,1e-8,1e-12]
        tab = []
        for scheme in schemes
            min_t = []
            steps = 0
            global mulc
            for i in 1:(reps)
                psi0 = setPsi0(t0)
                mulc  = 0
                line = @timed adaptive_time_stepper(psi0, t0, tend, dt, tol, scheme)
                min_t = push!(min_t, line[2])
                steps = line[1]
            end
            push!(tab, (steps, mulc, minimum(min_t), maximum(min_t), mean(min_t)))
        end
        @printf("%s %30s & \\#steps & \\#rejcs &  \\#mul! &       min time (s) &       max time (s) &      mean time (s) \\\\ \n","%","Scheme / $tol" )
        # println(tab)
        for i in 1:length(schemes)
            @printf("%s %30s & %7i & %7i & %7i & %18.12f & %18.12f & %18.12f \\\\ \n","%",schemes[i].name, tab[i][1][1],tab[i][1][2], tab[i][2], tab[i][3], tab[i][4], tab[i][5] )
        end
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        @printf("\\hline %30s & \\#steps &  time (s) \\\\ \\hline","Scheme / tol=$tol" )
        for i in 1:length(schemes)
            if testcase == "Rosen"
                @printf("\n %30s & %7i &  %6.4f \\\\",schemes[i].name, tab[i][1][1], tab[i][5] )
            else
                @printf("\n %30s & %7i &  %6.3f \\\\",schemes[i].name, tab[i][1][1], tab[i][5] )
            end

        end
        @printf(" \\hline \n")
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    end
end

#the first run compiles the benchmark
println("--------------------------------------Compiling $testcase--------------------------------------")
bm(t0, tend, dt, 1, testcase)
#now do $reps iterations
println("--------------------------------------Start BM $testcase --------------------------------------")
bm(t0,tend,dt,reps,testcase)
