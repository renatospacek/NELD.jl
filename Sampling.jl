module Sampling

using Polynomials
using Setup
using Integrators
using Compute
using KR
using Plots

export weakSampling!,
       weakOrder,
       clearSampling!,
       orderFit

# ==============================================================================
#                         Error sampling functions
# ==============================================================================
function weakSampling!(X::Array{Particle, 1},
                         Y::GenKR,
                         Z::Clist,
                         λ::Params,
                         C::Conv,
                         Niter::Int64,
                         fxn::Function)

SOILEA!(X[1], Y, Z, λ, λ.dt)
PBC!(X[1], Y, λ)
fill!(X[1].dW1, 0.0)
fill!(X[1].dW2, 0.0)

    for k in 2:C.p
        R = 2^(k-1)
        if mod(Niter, R) == 0
            fxn(X[k], Y, Z, λ, R*λ.dt)
            PBC!(X[k], Y, λ)
            fill!(X[k].dW3, 0.0)
            fill!(X[k].dW4, 0.0)
        end
    end
end

function weakOrder(X::Array{Particle,1}, Y::GenKR, C::Conv)
    for i in 1:C.p
        addKE(X[i])
        C.obs[i] += X[i].ke
    end
end

function clearSampling!(X::Array{Particle,1}, C::Conv)
    for k in 1:C.p
        fill!(X[k].q, 0.0)
        fill!(X[k].p, 0.0)
        fill!(X[k].qtmp, 0.0)
        fill!(X[k].prel, 0.0)
        fill!(X[k].dW, 0.0)
        X[k].ke = 0.0
        X[k].pe = 0.0
    end
end

function orderFit(C::Conv)
    C.obs /= C.M

    for i in 2:C.p
        C.err[i-1] = abs(C.obs[1] - C.obs[i])
    end

    xi = 2
    xf = C.p - 1
    yi = xi - 1
    yf = xf - 1

    x = C.h[xi:xf]
    y = C.err[yi:yf]

    p1 = polyfit(log.(x), log.(y), 1)
    c = coeffs(p1)

    #println("Order of weak convergence = $(c[2]) \n")

    return x, y
end



end
