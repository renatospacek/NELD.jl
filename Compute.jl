module Compute

using Setup
using LinearAlgebra

export addKE,
       fLJ,
       peLJ


# ==============================================================================
#                               Computing KE
# ==============================================================================
function addKE(X::Particle)
    X.ke = dot(X.prel,X.prel)/2/X.nPart
end

# ==============================================================================
#                          LJ force and potential
# ==============================================================================
function fLJ(rr::Float64, λ::Params)
    if rr > λ.rc
        return 0.0
    else
        r1 = λ.σ/rr
        r2 = r1^2

        return 48*λ.ε*r2^3*(r2^3 - 1/2)*r1
    end
end

function peLJ(rr::Float64, λ::Params)
    if rr > λ.rc
        return 0.0
    else
        r1 = λ.σ/rr
        r2 = r1^2

        return 4*λ.ε*r2^3*(r2^3 - 1) + 1
    end
end

end
