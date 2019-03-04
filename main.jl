#!/usr/local/bin/julia

include("input.jl")
include("initialize.jl")
include("integrator.jl")
include("histogram.jl")

# ==============================================================================
#                           Execute simulation
# ==============================================================================


function main()
    fxn = getfield(Main, integ)
    
    X, Y = initialize_box()
    
    X = compute_force(X, Y)
    kinetic!(X)
    
    L0 = zeros(Float64, 3)
    L0[1] = norm(Y.L[1,:])
    L0[2] = norm(Y.L[2,:])
    L0[3] = norm(Y.L[3,:])
    
    println("Size of simulation box: $(L0)")

    dW = sqrt(dt)*randn(dim, nPart*N)
    Nsteps = trunc(Int, N/R)

    for i in 1:Nsteps
        X = fxn(R, dW, X, i, Y)
        kinetic!(X)
        #X = pressure(X)
        
        if mod(i, outFreq) == 0
            println("Step $i | Max. F = $(maximum(X.f)) | Max. p = $(maximum(X.p)) | Max. prel = $(maximum(X.prel)) | PE = $(X.pe) | KE = $(X.ke) | T = $(X.temp) | Tot. E $(X.ke + X.pe)") 
           
        end
    end
    
    histogram(X.sp)
end




# ==============================================================================
@time main()








