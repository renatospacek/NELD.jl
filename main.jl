#!/usr/local/bin/julia

include("input.jl")
include("initialize.jl")
include("integrator.jl")
include("histogram.jl")

#using PyCall
#pygui(:qt)
#using PyPlot

#using PyPlot
#using Profile
#using ProfileView

# ==============================================================================
#                           Execute simulation
# ==============================================================================


function main()
   # fxn = getfield(Main, integ)
    
    X, Y = initialize_box()
    compute_force!(X, Y)
    kinetic!(X)
    

    dW = sqrt(dt)*randn(dim, nPart*N)
    Nsteps = trunc(Int, N/R)

    for i in 1:Nsteps
        EM!(R, X,dW, i, Y)
        kinetic!(X)
        #X = pressure(X)
        
        if mod(i, outFreq) == 0
        println("step $i")
        #    println("Step $i | Max. F = $(maximum(X.f)) | Max. p = $(maximum(X.p)) | T = $(X.temp)")        
        end
    end
    
    histogram(X.sp)
end




#function benchmark()
#    # Any setup code goes here.

#    # Run once, to force compilation.
#    println("======================= First run:")
#   # srand(666)
#    @time main()

#    # Run a second time, with profiling.
#    println("\n\n======================= Second run:")
#   # srand(666)
#    Profile.init(delay=0.01)
#    Profile.clear()
#    Profile.clear_malloc_data()
#    @profile @time main()

#    # Write profile results to profile.bin.
#    r = Profile.retrieve()
#    ProfileView.view()
#    #f = open("profile.bin", "w")
#    #Profile.print(f)
#   
#end














# ==============================================================================
@time main()
#Profile.init(delay=0.01)
#Profile.clear()
#Profile.clear_malloc_data()
#@profile @time main()
#f = open("profile2.bin", "w")
#Profile.print(f)
#ProfileView.view()
#benchmark()








