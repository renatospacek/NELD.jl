#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
function step_deform(Y::genkr)
    Y.θ += Y.δ*dt
    Y.θ -= round.(Y.θ)
    Y.εt = Y.ω*Y.θ
    expA = exp(Diagonal(Y.εt))
    Y.L = S*expA*Y.Vinv
    Y.Linv = Y.V*exp(Diagonal(-Y.εt))*Sinv
    
    return Y
end


# ==============================================================================
function PBC(X::particle, Y::genkr)
    qtmp = Array{Float64, 2}(undef, 3, nPart)
    qtmp = Y.Linv*X.q

	for i in 1:dim, j in 1:nPart
	    qtmp[i,j] -= round(qtmp[i,j])
	end
	
	X.q = Y.L*qtmp
	X.p = X.prel + A*X.q
	
	return X
end

# ==============================================================================
function lengthBC(dr::Array{Float64}, Y::genkr)
    ddr = zeros(Float64, dim)
    tdcoord = zeros(Float64, dim)
    dxin = copy(dr)

    for i in 1:3
	    ddr[i] = Y.Linv[i,1]*dxin[1] + Y.Linv[i,2]*dxin[2] + Y.Linv[i,3]*dxin[3]
	end
	
	for i in 1:3
	    tdcoord[i] = ddr[i] - round(ddr[i])
	end

    for i in 1:3
    	ddr[i] = Y.L[i,1]*tdcoord[1] + Y.L[i,2]*tdcoord[2] + Y.L[i,3]*tdcoord[3]
	end
	
	dist = norm(ddr)
	ddr = ddr/dist
	
	return ddr, dist
end




