#!/usr/local/bin/julia

include("input.jl")
include("initialize.jl")
include("integrator.jl")
include("histogram.jl")

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



# ==============================================================================

#function benchmark()
#
#    println("======================= First run:")
#    @time main()
#
#    println("\n\n======================= Second run:")
# 
#    Profile.init(delay=0.01)
#    Profile.clear()
#    Profile.clear_malloc_data()
#    @profile @time main()
#
#  
#    r = Profile.retrieve()
#    ProfileView.view()
#    #f = open("profile.bin", "w")
#    #Profile.print(f)
#   
#end

# ==============================================================================
@time main()
#benchmark








