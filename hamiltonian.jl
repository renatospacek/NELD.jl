#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#		            Lennard-Jones Potential and Force
# ==============================================================================
function fLJ(rr::Float64)::Float64
    if rr < rc
        r1 = σ/rr
        r2 = r1^2
    
        ff = 48*ε*r2^3*(r2^3 - 1/2)*r1
    else
        ff = 0.0
    end
    
	return ff
end

# ==============================================================================
function peLJ(rr::Float64)::Float64
    if rr < rc
        r1 = σ/rr
        r2 = r1^2
    
        pe = 4*ε*r2^3*(r2^3 - 1) + 1
    else
        pe = 0.0
    end
    
	return pe
end

# ==============================================================================
#                       Kinetic energy computation
# ==============================================================================
function kinetic!(X::particle{Float64})

    p = Array{Float64}(undef,3,nPart)
    
    for i in 1:3, j in 1:nPart
        @inbounds p[i,j] = X.prel[i,j]*X.prel[i,j]
    end

    p = vec(sum(p, dims=1))
    X.ke = sum(p)/2
    X.temp = X.ke*tfac
    storesp!(X.sp, p)
end


# ==============================================================================
function storesp!(sp::Array{Float64,1}, p::Array{Float64,1})

    for i in eachindex(p)
        p[i] = sqrt(p[i])
    end
    
    append!(sp,p)
end

# ==============================================================================
#                          Pressure Tensor
# ==============================================================================
function pressure!(X::particle{Float64})

	X.press = Array{Float64}(undef, 3, 3)
	
	for j in 1:dim, k in 1:dim, i in 1:nPart
	     @inbounds X.press[j,k] += X.prel[j,i]*X.prel[k,i]
	end
	
	@. X.press *= 1/vol

end






