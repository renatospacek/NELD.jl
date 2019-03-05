#!/usr/local/bin/julia

include("input.jl")
include("hamiltonian.jl")
include("kr.jl")

# ==============================================================================
#	                         Initialize list
# ==============================================================================

function build_list!(list::Array{Int64,1}, 
                     head::Array{Int64,1},
                     M::Array{Int64,1},
                     X::particle{Float64},
                     Y::genkr{Float64})

    mc = Array{Int64}(undef,3)
    da = Array{Float64}(undef,3)
    
    celn!(M, Y)
    
    for i in 1:nPart
        
        vec_index!(mc, da, Y, X, M, i)
        
        c = 1 + mc[1]*M[2]*M[3] + mc[2]*M[3] + mc[3]
        
        list[i] = head[c]
        head[c] = i
    end
end

# ==============================================================================
function clear_list!(list::Array{Int64,1},
                     head::Array{Int64,1})

    fill!(list, 0)
    fill!(head, 0)

end

# ==============================================================================
function compute_force!(X::particle{Float64},
                        Y::genkr{Float64})
                        
    fill!(X.f, 0.0)
    X.pe = 0.0
    
    if nopot == 1
        return X
    end
    
    list = Array{Int64,1}(undef, nPart) 
    head = Array{Int64,1}(undef, Mmax^3) 
    M = Array{Int64, 1}(undef, 3)
    
    clear_list!(list, head)
    build_list!(list, head, M, X, Y)
    
    mc = Array{Int64}(undef,3)
    ci = Array{Int64}(undef,3)
    da = Array{Float64}(undef,3)
    
    for i in 1:nPart
        vec_index!(mc, da, Y, X, M, i)
        
        for ci[1] in (mc[1]-1):(mc[1]+1),
            ci[2] in (mc[2]-1):(mc[2]+1),
            ci[3] in (mc[3]-1):(mc[3]+1)

            adjc = cell_index(ci, M)
            j = head[adjc]
             
            while j != 0
                if i<j   
                     update_force!(X, Y, i, j)
                end
                
                @inbounds j = list[j]
            end
        end
    end
end

# ==============================================================================
function update_force!(X::particle{Float64}, 
                       Y::genkr{Float64},
                       i::Int64,
                       j::Int64)


    dr = Array{Float64,1}(undef,3)
    r = Array{Float64,1}(undef,3)
    
    for l in 1:3
        dr[l] = X.q[l,i] - X.q[l,j]
    end
    
    dist = lengthBC(dr, r, Y)   
    ff = fLJ(dist)
    X.pe += peLJ(dist)

    for l in 1:3    
        X.f[l,i] += ff*r[l]
        X.f[l,j] -= ff*r[l]
    end        
end

# ==============================================================================
function vec_index!(mc::Array{Int64,1}, 
                    da::Array{Float64,1},
                    Y::genkr{Float64},
                    X::particle{Float64},
                    M::Array{Int64,1},
                    i::Int64)
    
    for j in 1:3
       @inbounds da[j] = Y.Linv[j,1]*X.q[1,i] + Y.Linv[j,2]*X.q[2,i] + Y.Linv[j,3]*X.q[3,i] + 0.5
       @inbounds mc[j] = floor(Int64, da[j]*M[j])
        
        if mc[j] < 0 || mc[j] > M[j]
            error("Bad index for particle $(i).")
        end
    end

end



# ==============================================================================
function celn!(M::Array{Int64,1}, Y::genkr{Float64})

    nL = Array{Float64}(undef, 3)
    
    cL = cofmat(Y.L)
    for i in 1:3
        nL[i] = norm(cL[:,i])*vol/rc

        M[i] = min(nL[i], Mmax)
        
        if M[i] < 3
            error("Less than 3 boxes in direction $(i), cell list invalid.")
        end
    end
end


# ==============================================================================
function cofmat(A::Array{Float64,2})::Array{Float64,2}

    cofA = Array{Float64}(undef, size(A))
    
    @. cofA = det(A)*inv(A)
    
    return cofA'
end


# ==============================================================================
function cell_index(ci::Array{Int,1}, M::Array{Int,1})::Int64

    c1 = mod(ci[3], M[3])
    c2 = mod(ci[2], M[2])*M[3]
    c3 = mod(ci[1], M[1])*M[2]*M[3]
    
    c = 1 + c1 + c2 + c3
    
    return c
end


