#!/usr/local/bin/julia

using Printf
using LinearAlgebra

# ==============================================================================
#                           Simulation parameters
# ==============================================================================

const t 	        = 10.0                               
const N 	        = 10^5                    
const dt            = t/N                                  
const dim	        = 3                                   
const ndim          = 5                                  
const nPart         = ndim^dim                                                
const γ 	        = 1.0                                  
const β 	        = 1.0                                  
const σ	            = 1.0                                 
const ε	            = 1.0   
const nopot         = 0                           
const rc	        = 2^(1/6)*σ                           
const R             = 1        
const outFreq       = 100
const η             = 0.0
const a             = 1.0
A                   = η*[1.0 0 0; 0 -1/2 0; 0 0 -1/2]
S                   = [1.0 0 0; 0 1 0; 0 0 1]*Diagonal([a*ndim; a*ndim; a*ndim])
Sinv                = inv(S) 
integrator          = "EM"


const vol = a^3*nPart
const Mmax = floor(Int, a/rc*ndim)

integ = Symbol(integrator)

# ==============================================================================
#                         Base/Type definitions
# ==============================================================================

Base.show(io::IO, f::Float64) = @printf(io, "%1.2f", f)

mutable struct particle{T<:AbstractFloat}
	q::Array{T}
	qtmp::Array{T}
	p::Array{T}
	prel::Array{T}
	f::Array{T}
	sp::Array{T}
	pe::Float64
	ke::Float64
	temp::T
	press::Array{T}
	function particle(nPart::Int64, β::T) where T<:AbstractFloat
	    q = zeros(Float64, 3,nPart)
	    p = zeros(Float64, 3,nPart)
	    qtmp = zeros(Float64, 3,nPart)
	    prel = zeros(Float64, 3,nPart)
	    f = zeros(Float64, 3, nPart)
	    sp = Array{T,1}(undef,nPart)
	    pe::T = 0.0
	    ke::T = 0.0
	    temp::T = 1/β
	    press = Array{T,2}(undef,3,3)   
	    new{T}(q,qtmp,p,prel,f,sp,pe,ke,temp,press)
    end
end

mutable struct genkr{T<:AbstractFloat}
    θ::Array{T}
    εt::Array{T}
    L::Array{T}
    Linv::Array{T}
    V::Array{T}
    Vinv::Array{T}
    ω::Array{T}
    δ::Array{T}
    function genkr(θ::Array{T}, S::Array{T}, Sinv::Array{T}) where T<:AbstractFloat
#        V = zeros(Float64, 3,3)
#        ω = zeros(Float64, 3,2)
#        
#        tmp1 = 0.73697622909957805;
#		tmp2 = 0.59100904850610347;
#		tmp3 = 0.32798527760568191;
#		V[1,1] = tmp2; V[1,2] = tmp3; V[1, 3] =-tmp1
#		V[2,1] =-tmp1; V[2,2] = tmp2; V[2, 3] =-tmp3
#		V[3,1] = tmp3; V[3,2] = tmp1; V[3, 3] = tmp2
#        
#        ω[1,1]=  1.619173832; ω[1,2]= -1.177725212
#		ω[2,1]= -0.441448621; ω[2,2]=  1.619173832 
#		ω[3,1]= -1.177725212; ω[3,2]= -0.441448621
#        
#        Vinv = V'
         
        M1 = [1.0 1 1; 1 2 2; 1 2 3]
        M2 = [2.0 -2 1; -2 3 -1; 1 -1 1]
        
        V = eigvecs(M1); Vinv = V'
        ω1 = log.(diag(Vinv*M1*V))

        V2 = eigvecs(M2); V2inv = V2'
        ω2 = log.(diag(V2inv*M2*V2))
        
        ω = [ω1 ω2]
        δ = ω\diag(A)
        εt = ω*θ
        expA = exp(Diagonal(εt))
        L = S*expA*Vinv
        Linv = inv(L)
        new{T}(θ,εt,L,Linv,V,Vinv,ω,δ)
    end
end





