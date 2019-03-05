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


## ==============================================================================
#function split(R::Int64, dW::Array{Float64,2}, X::particle, i::Int, Y::genkr) 
#    h = dt*R
#	Winc = _addnoise(R, dW, i)
#	
#	X.p += 1/2*X.f*h
#	X.q += X.p*h
#	X.prel = X.p - A*X.q
#	
#	Y = step_deform(Y)
#    X = PBC(X, Y)
#    X = compute_force(X, Y)
#    
#    rf = randn(3,nPart)
#    
#    decay = exp(-γ*h)
#    diffuse = sqrt((1-decay^2)/β)
#    
#    X.p += 1/2*X.f*h
#	X.p += A*X.p*h
#	
#	X.p = X.p*decay + A*X.q*(1-decay) + rf*diffuse
#	X.prel = X.p - A*X.q
#	
#	kinetic!(X)
#	
#	return X
#end

## ==============================================================================
#function ABAPO(R::Int64, dW::Array{Float64,2}, X::particle, i::Int, Y::genkr) 
#    h = dt*R
#	Winc = _addnoise(R, dW, i)
#	
#	Y = step_deform(Y)
#    X = PBC(X, Y)
#    #X = compute_force(X, Y)
#    
#	X.q += X.p*h
#	X.prel = X.p - A*X.q
#	
#	Y = step_deform(Y)
#    X = PBC(X, Y)
#    #X = compute_force(X, Y)
#    
#    rf = randn(3,nPart)
#    
#    decay = exp(-γ*h)
#    diffuse = sqrt((1-decay^2)/β)
#    
#    for i in 1:nPart
#	    X.p[:,i] = exp(A*dt)*X.p[:,i]
#	end
#	
#	X.p = γ*X.p + A*X.q*(1-γ) + rf*diffuse
#	X.prel = X.p - A*X.q
#	
#	return X
#end






## ==============================================================================
#function _addnoise(R::Int, dW::Array{Float64,2}, j::Int)
#    Winc = Array{Float64, 2}(undef, dim, nPart)
#    
#    for l in 1:dim, i in 1:nPart
#        Winc[l,i] = sum(dW[l,(nPart*R*(j-1)+R*(i-1)+1):R*(nPart*(j-1)+i)])
#    end
#    
#    return Winc
#end

# k==============================================================================
#                           Symplectic Euler A
# ==============================================================================

#function SEA(R::Int64, dW::Array{Float64,2}, X::particle, i::Int) 
#    h = dt*R
#	Winc = _addnoise(R, dW, i)
#    
#    # Apply PBC - GenKR(qi, pi)
#	X.q +=  h*X.p
#	# Apply PBC - GenKR(q(i+1),pi)
#	X.f = (LJ(r) - γ*(X.p - A*X.q) + A*X.q)*h
#	X.p +=  h*X.f + sqrt(2*γ/β)*Winc
#	
#	return X, pe
#end



# ==============================================================================
#                     Compute additive noise array
# ==============================================================================

function _addnoise(R::Int, dW::Array{Float64}, j::Int)
   # Winc = Array{Float64}(undef,3,nPart)
    
    #for i in 1:nPart
        #Winc[:,i] = sum(dW[:,(nPart*R*(j-1)+R*(i-1)+1):R*(nPart*(j-1)+i)], dims=2)
        #Winc[:,i] .= dW[:, nPart*(j-1)+i]
        
   # end
    
    return view(dW, :, (nPart*(j-1)+1):j*nPart)   #Winc
end















