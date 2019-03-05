#!/usr/local/bin/julia

using Printf
using LinearAlgebra

# ==============================================================================
#                           Simulation parameters
# ==============================================================================
const t 	        = 10.0                               
const N 	        = 5*10^4                   
const dt            = t/N                                  
const dim	        = 3                                   
const ndim          = 10                                  
const nPart         = ndim^dim                                                
const γ 	        = 1.0                                  
const β 	        = 1.0                                  
const σ	            = 1.0                                 
const ε	            = 1.0   
const nopot         = 1                           
const rc	        = 2^(1/6)*σ                           
const R             = 1        
const outFreq       = 1000
const η             = 0.5
const a             = 1.0
const Σ             = sqrt(2*γ/β)
 S             = a*ndim*Array{Float64}(I, 3, 3)
 Sinv          = inv(S) 
 A             = η*[1.0 0 0; 0 -1/2 0; 0 0 -1/2]
#integrator          = "EM"

const tfac = 2/3/nPart
const vol = a^3*nPart
const Mmax = floor(Int, a/rc*ndim)

    
    
#integ = Symbol(integrator)

# ==============================================================================
#                         Base/Type definitions
# ==============================================================================

Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)

function Base.exp(A::Array{T,1}) where T<:Real
    return exp(Diagonal(A))
end

#eye(T::Type, N::Int) = Array{T,2}(I, N, N)
#eye(N::Int) = Array{Float64,2}(I, N, N)















mutable struct particle{T<:AbstractFloat}
	q::Array{T,2}
	p::Array{T,2}
	prel::Array{T,2}
	f::Array{T,2}
	sp::Array{T,1}
	pe::Float64
	ke::Float64
	temp::Float64
	press::Array{T,2}
	function particle(nPart::Int64, β::T) where T<:AbstractFloat
	    q = zeros(Float64, 3,nPart)
	    p = zeros(Float64, 3,nPart)
	    prel = zeros(Float64, 3,nPart)
	    f = zeros(Float64, 3, nPart)
	    sp = Array{T,1}(undef,nPart)
	    pe::T = 0.0
	    ke::T = 0.0
	    temp::T = 1/β
	    press = Array{T,2}(undef,3,3)   
	    
	    new{T}(q,p,prel,f,sp,pe,ke,temp,press)
    end
end

mutable struct genkr{T<:AbstractFloat}
    θ::Array{T,1}
    εt::Array{T,1}
    L::Array{T,2}
    Linv::Array{T,2}
    V::Array{T,2}
    Vinv::Array{T,2}
    ω::Array{T,2}
    δ::Array{T,1}
    function genkr(θ::Array{T,1}, S::Array{T,2}, Sinv::Array{T,2}) where T<:AbstractFloat
        M1 = [1.0 1 1; 1 2 2; 1 2 3]
        #M2 = [2.0 -2 1; -2 3 -1; 1 -1 1]
        
        V = eigvecs(M1); Vinv = V'
        ω1 = log.(diag(Vinv*M1*V))

        ω2 = [ω1[2]; ω1[3]; ω1[1]]
        
        ω = [ω1 ω2]
        δ = ω\diag(A)
        εt = ω*θ
        L = S*exp(εt)*Vinv
        Linv = inv(L)
        
        new{T}(θ,εt,L,Linv,V,Vinv,ω,δ)
    end
end





