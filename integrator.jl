#!/usr/local/bin/julia

include("input.jl")
include("cell_lists.jl")

# ==============================================================================
#                             Euler-Maruyama
# ==============================================================================

function EM!(R::Int64, 
             X::particle{Float64},
             dW::Array{Float64,2},
             i::Int,
             Y::genkr{Float64})
             
             
	qtmp = Array{Float64}(undef, 3, nPart)
	h = dt*R
	Winc = _addnoise(R, dW, i)
	
	step_deform!(Y)
        PBC!(X, Y)
	compute_force!(X, Y)
    
	qtmp = X.q + h*X.p
	F = (X.f - γ*(X.p - A*X.q) + A*X.p)*h
	X.p +=  F + Σ*Winc
	X.q = qtmp
	
	X.prel = X.p - A*X.q
end

# ==============================================================================
#                     Compute additive noise array
# ==============================================================================

function _addnoise(R::Int, dW::Array{Float64,2}, j::Int)
    
    return view(dW, :, (nPart*(j-1)+1):j*nPart)
end






