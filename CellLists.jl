module CellLists

using LinearAlgebra
using Setup
using KR
using Compute

export computeForce

# ==============================================================================
#	                    Construct cell list
# ==============================================================================
function build_list!(X::Particle, Y::GenKR, Z::Clist, λ::Params)
    celn!(Y, Z, λ)

    for i in 1:λ.nPart
        vec_index!(X, Y, Z, i)
        Z.list[i] = Z.head[Z.c]
        Z.head[Z.c] = i
    end
end

# ==============================================================================
function clear_list!(Z::Clist)
    fill!(Z.list, 0)
    fill!(Z.head, 0)
end

# ==============================================================================
function computeForce(X::Particle,
                        Y::GenKR,
                        Z::Clist,
                        λ::Params)

    fill!(X.f, 0.0)
    X.pe = 0.0

    if λ.LJ == false
        return X
    end

    clear_list!(Z)
    build_list!(X, Y, Z, λ)
     for i in 1:λ.nPart
        vec_index!(X, Y, Z, i)
         for ci1 in (Z.mc[1]-1):(Z.mc[1]+1),
             ci2 in (Z.mc[2]-1):(Z.mc[2]+1),
             ci3 in (Z.mc[3]-1):(Z.mc[3]+1)

            Z.c = cell_index(ci1,ci2,ci3, Z.M)
            j = Z.head[Z.c]

            while j != 0
                if i<j
                     updateForce!(X, Y, Z, λ, i, j)
                end
                @inbounds j = Z.list[j]
            end
        end
    end
end

# ==============================================================================
function updateForce!(X::Particle,
                       Y::GenKR,
                       Z::Clist,
                       λ::Params,
                       i::Int64,
                       j::Int64)

    Z.dr[1] = X.q[1,i] - X.q[1,j]
    Z.dr[2] = X.q[2,i] - X.q[2,j]
    Z.dr[3] = X.q[3,i] - X.q[3,j]

    dist = lengthBC(Y, Z)
    ff = fLJ(dist, λ)
    X.pe += peLJ(dist, λ)

    @inbounds for l in 1:3
        X.f[l,i] += ff*Z.r[l]
        X.f[l,j] -= ff*Z.r[l]
    end
end

# ==============================================================================
function vec_index!(X::Particle,
                    Y::GenKR,
                    Z::Clist,
                    i::Int64)


    Z.da .= Y.Linv*X.q[:,i] .+ 0.5*ones(3)
    Z.mc .= floor.(Int64, Z.da.*Z.M)
    for i in 1:3
        Z.mc[i] = mod(Z.mc[i], Z.M[i])
    end

    for j in 1:3
        if Z.mc[j] < 0 || Z.mc[j] > Z.M[j]
            error("Bad index for particle $(i). \n Cell mc[$j] = $(Z.mc[j]).")
        end
    end
    Z.c = cell_index(Z.mc[1], Z.mc[2], Z.mc[3], Z.M)
    if Z.c > (Z.Mmax)^3
        error("Nonexistence cell index $(Z.c)")
    end
end

# ==============================================================================
function celn!(Y::GenKR, Z::Clist, λ::Params)
    cL = cofmat(Y.L)

    @inbounds for i in 1:3
        Z.nL[i] = λ.fac./norm(cL[:,i])

        Z.M[i] = floor(Int64, min(Z.nL[i], λ.Mmax))
        if Z.M[i] < 3 || Z.M[i] > λ.Mmax
            error("Less than 3 boxes in direction $(i), cell list invalid.")
        end
    end
end
# ==============================================================================
function cell_index(ci1::Int64, ci2::Int64, ci3::Int64, M::Array{Int,1})
    c1::Int64 = wrap(ci3, M[3])
    c2::Int64 = wrap(ci2, M[2])*M[3]
    c3::Int64 = wrap(ci1, M[1])*M[2]*M[3]

    return 1 + c1 + c2 + c3
end
# ==============================================================================
function cofmat(A::Array{Float64,2})
    cofA = Array{Float64}(undef, size(A))
   # @. cofA = det(A)*inv(A)
    for i in 0:2, j in 0:2
        cofA[i+1,j+1] = A[(i+1)%3+1, (j+1)%3+1]*A[(i+2)%3+1, (j+2)%3+1]
    end
    return cofA
end

end
