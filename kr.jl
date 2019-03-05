#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
#                         KR Boundary Conditions
# ==============================================================================
function step_deform!(Y::genkr{Float64})

    Y.θ[1] += Y.δ[1]*dt
    Y.θ[2] += Y.δ[2]*dt
    
    Y.θ[1] -= round(Y.θ[1])
    Y.θ[2] -= round(Y.θ[2])
    
    Y.εt = Y.ω*Y.θ
    
    Y.L = S*exp(Y.εt)*Y.Vinv
    Y.Linv = inv(Y.L) #Y.V*exp(-Y.εt)/a/ndim
end


# ==============================================================================
function PBC!(X::particle{Float64}, Y::genkr{Float64})

    qtmp = Array{Float64, 2}(undef, 3, nPart)
    qtmp = Y.Linv*X.q

	for i in 1:dim, j in 1:nPart
	    qtmp[i,j] -= round(qtmp[i,j])
	end
	
	X.q = Y.L*qtmp
	X.p = X.prel + A*X.q
	
end

# ==============================================================================
function lengthBC(dr::Array{Float64,1}, r::Array{Float64,1}, Y::genkr{Float64})

    tdcoord = Array{Float64}(undef, dim)

    for i in 1:3
	    r[i] = Y.Linv[i,1]*dr[1] + Y.Linv[i,2]*dr[2] + Y.Linv[i,3]*dr[3]
	end
	
	for i in 1:3
	    tdcoord[i] = r[i] - round(r[i])
	end

    for i in 1:3
    	r[i] = Y.L[i,1]*tdcoord[1] + Y.L[i,2]*tdcoord[2] + Y.L[i,3]*tdcoord[3]
	end
	
	dist::Float64 = norm(r)
	normalize!(r)
	
	return dist
end




