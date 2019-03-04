#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#		            Lennard-Jones Potential and Force
# ==============================================================================
function LJ(rr::Float64)
    
    if rr < rc
        r1 = σ/rr
        r2 = r1^2
    
        ff = 48*ε*r2^3*(r2^3 - 1/2)*r1
        pe = 4*ε*r2^3*(r2^3 - 1) + 1
    else
        ff = 0.0
        pe = 0.0
    end
    
	return ff, pe
end

# ==============================================================================
#                       Kinetic energy computation
# ==============================================================================
function kinetic!(X::particle)
    p = Array{Float64}(undef,3,nPart)
    
    for i in 1:3, j in 1:nPart
        @inbounds p[i,j] = X.prel[i,j]*X.prel[i,j]
    end
    
    p = vec(sum(p, dims=1))
    X.ke = sum(p)/2
    X.temp = X.ke*2/3/nPart
    append!(X.sp, sqrt.(p))
    
    return X
end


# ==============================================================================
#                          Pressure Tensor
# ==============================================================================
function pressure(X::particle)
	#X.press = zeros(Float64, dim, dim)
	
	for j in 1:dim, k in 1:dim, i in 1:nPart
	     X.press[j,k] += X.prel[j,i]*X.prel[k,i]
	end
	
	X.press *= 1/vol
	
	return X 
end






