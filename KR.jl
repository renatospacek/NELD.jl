module KR

using Setup
using LinearAlgebra

export stepDeform!,
       PBC!,
       lengthBC

# ==============================================================================
#                           Periodic Boundaries
# ==============================================================================
function stepDeform!(Y::GenKR, λ::Params, h::Float64)
    Y.θ[1] += Y.δ[1]*h
    Y.θ[2] += Y.δ[2]*h

    Y.θ[1] -= round(Y.θ[1])
    Y.θ[2] -= round(Y.θ[2])

    #Y.εt = Y.ω*Y.θ

    for i in 1:3
        Y.εt[i] = Y.ω[i]*Y.θ[1] + Y.ω[i]*Y.θ[2]
    end

    Y.L = λ.S*exp(Y.εt)*Y.Vinv
    Y.Linv = Y.V*exp(-Y.εt)*λ.Sinv
end

function PBC!(X::Particle, Y::GenKR, λ::Params)
    X.qtmp2 .= Y.Linv*X.q

	for i in 1:3, j in 1:λ.nPart
	    X.qtmp2[i,j] -= round(X.qtmp2[i,j])
	end

    X.q = Y.L*X.qtmp2
    X.p .= X.prel .+ λ.A*X.q
end

function lengthBC(Y::GenKR, Z::Clist)::Float64
    Z.r[1] = Y.Linv[1,1]*Z.dr[1] + Y.Linv[1,2]*Z.dr[2] + Y.Linv[1,3]*Z.dr[3]
    Z.r[2] = Y.Linv[2,1]*Z.dr[1] + Y.Linv[2,2]*Z.dr[2] + Y.Linv[2,3]*Z.dr[3]
    Z.r[3] = Y.Linv[3,1]*Z.dr[1] + Y.Linv[3,2]*Z.dr[2] + Y.Linv[3,3]*Z.dr[3]

    Z.tdcoord[1] = Z.r[1] - round(Z.r[1])
    Z.tdcoord[2] = Z.r[2] - round(Z.r[2])
    Z.tdcoord[3] = Z.r[3] - round(Z.r[3])

	Z.r[1] = Y.L[1,1]*Z.tdcoord[1] + Y.L[1,2]*Z.tdcoord[2] + Y.L[1,3]*Z.tdcoord[3]
	Z.r[2] = Y.L[2,1]*Z.tdcoord[1] + Y.L[2,2]*Z.tdcoord[2] + Y.L[2,3]*Z.tdcoord[3]
	Z.r[3] = Y.L[3,1]*Z.tdcoord[1] + Y.L[3,2]*Z.tdcoord[2] + Y.L[3,3]*Z.tdcoord[3]

    dist = norm(Z.r)

	Z.r = Z.r/dist

	return dist
end

function BCnorm(Y::GenKR, BCdist::Array{Float64, 2})

    dist = zeros(Float64, 3, size(BCdist,2))

    for i in 1:size(BCdist,2)
        dist[:, i] = lengthBC(Y, BCdist[:,i])
    end

    return norm(dist)
end

end
