#!/usr/local/bin/julia

import PyPlot

include("input.jl")
include("integrator.jl")

# ==============================================================================
#                         Histogram plot
# ==============================================================================
function histogram(p::Array{Float64})
    
    fxn1(x) = (x.^2).*exp.(-β*x.^2/2)
    Z = simps(fxn1, 0.0, 100.0, 10^4)
    
    x = range(0, stop=6, length=1000)
    y = 1/Z*fxn1(x)
    
	PyPlot.plt[:hist](p, bins=100, 
	                     range=(0,10),
	                     density=1, 
	                     ec="black",
	                     label="Simulation data")
    
    PyPlot.plot(x,y, color="red", 
	                 linewidth=2.0, 
	                 linestyle="-",
	                 label="Gibbs dist.")
	
    PyPlot.xlabel("||v||")
    PyPlot.ylabel("PDF")
    PyPlot.title("BEF | η = $η, β = $β, γ = $γ")
    PyPlot.legend()
    
	PyPlot.show()
end

# ==============================================================================
#                Simpson's rule to normalize the distribution
# ==============================================================================
function simps(f::Function, a::Float64, b::Float64, n::Int)::Float64
    h = (b-a)/n
    s = f(a) + f(b)
    s += 4*sum(f(collect(1:2:n)*h))
    s += 2*sum(f(collect(2:2:n-1)*h))
    
    return h/3*s
end










