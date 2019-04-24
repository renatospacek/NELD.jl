module Noise

using Setup
using Random

export noise!

# ==============================================================================
#                             Gaussian white noise
# ==============================================================================
function noise!(X::Array{Particle,1}, λ::Params, C::Conv, Niter::Int64, l::Int64)
    Random.seed!(Niter*l)
    #rd = rand(1:10^6)
    #Random.seed!(rd)
    xx = randn(Float64, 3, λ.nPart)
    yy = randn(Float64, 3, λ.nPart)

    X[1].dW1 .= xx.*λ.Σ
    X[1].dW2 .= (xx .+ yy.*(1/sqrt(3))).*0.5.*λ.Σ2
    X[1].dW3 .= xx.*λ.Σ
    X[1].dW4 .= (xx .+ yy.*(1/sqrt(3))).*0.5.*λ.Σ2

    for k in 2:C.p
        X[k].dW4 += X[1].dW4 .+ X[k].dW3.*λ.dt
        X[k].dW3 += X[1].dW3
    end
end




end
