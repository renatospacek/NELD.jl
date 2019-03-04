#!/usr/local/bin/julia

include("input.jl")

# ==============================================================================
function initialize_box()
    qlam = zeros(Float64, 3, 1)
    
    X = particle(nPart, β)

    Y = genkr([0.0; 0.0], S, Sinv)

    for l in 0:ndim-1, i in 0:ndim-1, j in 0:ndim-1
        if i*ndim+j < ndim^2
            qlam[1,1] = (0.5 + i-0.5*ndim)/ndim
            qlam[2,1] = (0.5 + j-0.5*ndim)/ndim
            qlam[3,1] = (0.5 + l-0.5*ndim)/ndim
            ind = floor(Int, (i*ndim+j)*ndim+l)+1
            qlam = Y.L*qlam
            X.q[1,ind] = qlam[1,1]
            X.q[2,ind] = qlam[2,1]
            X.q[3,ind] = qlam[3,1]
        end
    end
    
    X = _initialize_vel(X, Y)
    
    return X, Y
end


# ==============================================================================
function _initialize_vel(X::particle, Y::genkr)
    X.p = A*X.q
    
    for i in 1:3, j in 1:nPart
        X.q[i,j] += 0.05*randn()
        X.p[i,j] += sqrt(1/β)*randn()
    end
    
    X.prel = X.p - A*X.q
    α = sqrt(1/β*dim*nPart/dot(X.prel,X.prel))
    X.p -= X.prel*(1-α)
    X.prel = X.p - A*X.q
    X = PBC(X,Y)
    
    return X
end







