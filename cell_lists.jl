#!/usr/local/bin/julia

include("input.jl")
include("hamiltonian.jl")
include("kr.jl")

# ==============================================================================
#	                         Initialize list
# ==============================================================================

function initialize_list(X::particle, Y::genkr)
    list = zeros(Int, nPart) 
    head = zeros(Int, Mmax^3) 
#    da = zeros(Float64, 3)
#    mc = zeros(Int, 3) 
#    M = zeros(Int, 3)
    mc = Array{Int}(undef,3)
    M = Array{Int}(undef,3)
    da = Array{Float64}(undef,3)
    
    cL = cofmat(Y.L)
    nL = vol/rc*[norm(cL[:,1]); norm(cL[:,2]); norm(cL[:,3])]
    
    for i in 1:3
        M[i] = max(min(nL[i], Mmax), 1.0)
        
        if M[i] < 3
            error("Less than 3 boxes in direction $(i), cell list invalid.")
        end
    end
    
    Ncel = prod(M)
    
    for i in 1:nPart
        for j in 1:3
            da[j] = Y.Linv[j,1]*X.q[1,i] + Y.Linv[j,2]*X.q[2,i] + Y.Linv[j,3]*X.q[3,i] + 0.5
        end
        
        for a in 1:3
            mc[a] = floor(Int, da[a]*M[a])
            
            if mc[a] < 0 || mc[a] > M[a]
                error("Bad index for particle $(i).")
            end
        end
        
        c = 1 + mc[1]*M[2]*M[3] + mc[2]*M[3] + mc[3]
        
        list[i] = head[c]
        head[c] = i
    end
    
    return head, list, M
end


# ==============================================================================
function compute_force(X::particle, Y::genkr)
    X.f = zeros(Float64, dim, nPart)
    X.pe = 0.0
    
    if nopot == 1
        return X
    end
    
    head, list, M = initialize_list(X, Y)
#    mc = zeros(Int, 3) 
#    ci = zeros(Int, 3) 
#    da = zeros(Float64, 3) 
    
    mc = Array{Int}(undef,3)
    ci = Array{Int}(undef,3)
    da = Array{Float64}(undef,3)
    
    for i in 1:nPart
        for j in 1:3
            da[j] = Y.Linv[j,1]*X.q[1,i] + Y.Linv[j,2]*X.q[2,i] + Y.Linv[j,3]*X.q[3,i] + 0.5
        end
        
        for a in 1:3
            mc[a] = floor(Int, da[a]*M[a])
        end
        
        for ci[1] in (mc[1]-1):(mc[1]+1),
            ci[2] in (mc[2]-1):(mc[2]+1),
            ci[3] in (mc[3]-1):(mc[3]+1)

            adjc = cell_index(ci, M)
            j = head[adjc]
             
            while j != 0
                if i<j
                    dr = X.q[:,i] - X.q[:,j]
                    r, dist = lengthBC(dr, Y)   
                    ff, pe = LJ(dist) 
                    X.f[:,i] += ff*r
                    X.f[:,j] -= ff*r
                    X.pe += pe
                end
                j = list[j]
            end
        end
    end

    return X
end




# The function compute_force2 doesn't use cell lists to compute the force
#function compute_force2(X::particle, Y::genkr)
#    X.f = zeros(dim, nPart)
#    X.pe = 0.0
#    
#    if nopot == 1
#        return X
#    end
#    
#    for i in 1:nPart-1, j in (i+1):nPart
#            dr = X.q[:,i] - X.q[:,j]
#            r, dist = lengthBC(dr, Y)   
#            ff, pe = LJ(dist) 
#            X.f[:,i] += ff*r
#            X.f[:,j] -= ff*r
#            X.pe += pe
#    end
#    
#    return X
#end






# ==============================================================================
function cofmat(A::Array{Float64,2})
    cofA = det(A)*inv(A)
    cofA = cofA'
    
    return cofA
end


function cell_index(ci::Array{Int}, M::Array{Int})
    c1 = mod(ci[3], M[3])
    c2 = mod(ci[2], M[2])*M[3]
    c3 = mod(ci[1], M[1])*M[2]*M[3]
    c = 1 + c1 + c2 + c3
    
    return c
end


