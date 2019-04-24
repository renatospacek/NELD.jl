module Setup

using Random
using Parameters
using LinearAlgebra
using Statistics
using Polynomials
using Plots

export Params,
       Particle,
       GenKR,
       Conv,
       Clist,
       wrap

# ==============================================================================
#                           Struct definitions
# ==============================================================================

@with_kw struct Params
    t::Float64
    N::Int64
    dt::Float64 = t/N
    ndim::Int64
    γ::Float64 = 1.0
    β::Float64 = 1.0

    a::Float64
    S::Array{Float64,2} = a*ndim*Array{Float64,2}(I, 3, 3)
    A::Array{Float64,2}
    Γ::Array{Float64,2} = γ*Array{Float64,2}(I, 3, 3) - A
    Γ2::Array{Float64,2} = γ*Array{Float64,2}(I, 3, 3) - A

    Sinv::Array{Float64,2} = inv(S)
    Σ::Float64 = sqrt(2*dt*γ/β)
    Σ2::Float64 = sqrt(2*(dt^3)*γ/β)
    nPart::Int64 = ndim^3
    α::Float64 = sqrt(1/β*3*nPart)

    LJ::Bool = true
    ε::Float64 = 1.0
    σ::Float64 = 1.0
    rc::Float64 = 2^(1/6)*σ
    Mmax::Int64 = floor(Int64, a/rc*ndim)
    vol::Float64 = a^3*nPart
    fac::Float64 = vol/rc
    plot::Bool
    soile::Bool = false
end

@with_kw mutable struct Particle
    nPart::Int64
    q::Array{Float64,2} = zeros(Float64, 3, nPart)
    qtmp::Array{Float64,2} = zeros(Float64, 3, nPart)
    qtmp2::Array{Float64,2} = zeros(Float64, 3, nPart)
    p::Array{Float64,2} = zeros(Float64, 3, nPart)
    prel::Array{Float64,2} = zeros(Float64, 3, nPart)
    f::Array{Float64,2} = zeros(Float64, 3, nPart)
    C1::Array{Float64,2} = zeros(Float64, 3, 3)
    C2::Array{Float64,2} = zeros(Float64, 3, nPart)
    Ftmp::Array{Float64,2} = zeros(Float64, 3, nPart)
    dW::Array{Float64,2} = zeros(Float64, 3, nPart)
    dW1::Array{Float64,2} = zeros(Float64, 3, nPart)
    dW2::Array{Float64,2} = zeros(Float64, 3, nPart)
    dW3::Array{Float64,2} = zeros(Float64, 3, nPart)
    dW4::Array{Float64,2} = zeros(Float64, 3, nPart)
    ke::Float64 = 0.0
    pe::Float64 = 0.0
end

mutable struct GenKR
    θ::Array{Float64,1}
    εt::Array{Float64,1}
    L::Array{Float64,2}
    Linv::Array{Float64,2}
    ω::Array{Float64,2}
    δ::Array{Float64,1}
    V::Array{Float64,2}
    Vinv::Array{Float64,2}
    function GenKR(θ::Array{Float64,1}, A::Array{Float64,2}, S::Array{Float64,2})
        M1 = [1.0 1 1; 1 2 2; 1 2 3]
        V = eigvecs(M1)
        Vinv = inv(V)
        ω1 = log.(diag(Vinv*M1*V))
        ω2 = [ω1[2]; ω1[3]; ω1[1]]
        ω = [ω1 ω2]
        δ = ω\diag(A)
        εt = ω*θ
        L = S*exp(Diagonal(εt))*Vinv
        Linv = inv(L)
        new(θ,εt,L,Linv,ω,δ,V,Vinv)
    end
end

mutable struct Conv
    p::Int64
    M::Int64
    h::Array{Float64,1}
    obs::Array{Float64,1}
    err::Array{Float64,1}
    plot::Bool
    function Conv(p::Int64, M::Int64, λ)
        h = zeros(Float64, p)
        for i in 1:p
            h[i] = λ.dt*2^(i-1)
        end
        obs = zeros(Float64, p)
        err = zeros(Float64, p-1)
        plot = λ.plot
        new(p, M, h, obs, err, plot)
    end
end

@with_kw mutable struct Clist
    nPart::Int64
    Mmax::Int64
    c::Int64 = 0
    list::Array{Int64,1} = Array{Int64,1}(undef, nPart)
    head::Array{Int64,1} = Array{Int64,1}(undef, Mmax^3)
    M::Array{Int64,1} = Array{Int64, 1}(undef, 3)
    mc::Array{Int64,1} = Array{Int64, 1}(undef, 3)
    da::Array{Float64,1} = Array{Float64}(undef, 3)
    tdcoord::Array{Float64,1} = Array{Float64}(undef, 3)
    dr::Array{Float64,1} = Array{Float64}(undef, 3)
    r::Array{Float64,1} = Array{Float64}(undef, 3)
    nL::Array{Float64,1} = Array{Float64}(undef, 3)
end

# ==============================================================================
#                      Base/other function definitions
# ==============================================================================
function Base.exp(A::Array{T,1}) where T<:Real
    return exp(Diagonal(A))
end

function wrap(x::Int64, n::Int64)::Int64
    while x >= n; x -= n; end
    while x < 0; x += n; end

    return x
end

end
