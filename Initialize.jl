module Initialize

using Setup
using CellLists
using KR
using LinearAlgebra
export initialize

# ==============================================================================
#                      Initialization and resetting
# ==============================================================================
function initialize(λ::Params, C::Conv)
    Y = GenKR([0.0; 0.0], λ.A, λ.S)
    X = Array{Particle, 1}(undef, C.p)
    Z = Clist(nPart = λ.nPart, Mmax = λ.Mmax)

    for i in 1:C.p
        X[i] = Particle(nPart = λ.nPart)
    end
    initialize(X[1], Y, λ)
    for i in 2:C.p
        X[i].q .= X[1].q
        X[i].p .= X[1].p
    end

    for i in 1:C.p
        X[i].prel .= X[i].p - λ.A*X[i].q
        α = λ.α/sqrt(dot(X[i].prel,X[i].prel))
        X[i].p -= X[i].prel*(1 - α)
        X[i].prel .= X[i].p - λ.A*X[i].q
        PBC!(X[i], Y, λ)
        computeForce(X[i], Y, Z, λ)
    end

    return X, Y, Z
end

function initialize(X::Particle, Y::GenKR, λ::Params)
    rg = 0:(λ.ndim-1)

    ll = 1
    for l in rg, i in rg, j in rg
        X.q[3,ll] = (0.5 + i-0.5*λ.ndim)/λ.ndim
        X.q[2,ll] = (0.5 + j-0.5*λ.ndim)/λ.ndim
        X.q[1,ll] = (0.5 + l-0.5*λ.ndim)/λ.ndim
        ll += 1
    end

    X.q .= Y.L*X.q
    X.q += 0.05*randn(3,λ.nPart)
    X.p .= λ.A*X.q
    X.p += sqrt(1/λ.β)*randn(3,λ.nPart)

end

end
