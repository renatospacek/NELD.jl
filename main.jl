using WeakNELD
using Parameters
using DelimitedFiles
using Random
using LinearAlgebra
using Plots
using Polynomials

# ==============================================================================
#                             Main function
#
# TO DO:
# - try out different boundary condition implementations (regular and corrected)
# ==============================================================================
function main()
    λ =  Params(t = 1.0/16,
                N = 2^10,
                ndim = 5,
                a = 1.25,
                LJ = false,
                plot = false,
                #A = [0.0 0 0; 0 0 0; 0 0 0])
                A = 0.1*[1.0 0 0; 0 -1/2 0; 0 0 -1/2])

    outFreq = 2^10
    C = Conv(10, 8, λ) # p, M
    Rmax = 2^(C.p-1)
    Nmin = Int(λ.N/Rmax)
    obs = zeros(Float64, Nmin, C.p)

    for l in 1:C.M
        println("M = $l of $(C.M) ==========")
        X, Y, Z = initialize(λ, C)

        for Niter in 1:λ.N
            noise!(X, λ, C, Niter, l)
            weakSampling!(X, Y, Z, λ, C, Niter, EM!)
            stepDeform!(Y, λ, λ.dt)

            if mod(Niter, outFreq) == 0
                println("$(Int(Niter/outFreq))/$(Int(λ.N/outFreq))")
            end

            if mod(Niter, Rmax) == 0
                for ii in 1:C.p
                    addKE(X[ii])
                    obs[Int(Niter/Rmax), ii] += X[ii].ke
                end
            end
        end

        weakOrder(X, Y, C)
        clearSampling!(X, C)
    end
    obs /= C.M
    Nvec = collect(Rmax:Rmax:λ.N)
    io = open("EMtest.txt","w")
    writedlm(io, [Nvec obs], " ")
    close(io)
    x, y = orderFit(C)

    #println("\n h = $(C.h)")
    #println("err = $(C.err)")
    #println("ke = $(C.obs)")

    if C.plot == true
        showPlot(x, y)
    end

    return x, y
end

## ==============================================================================

function showPlot(x::Array{Float64,1},
                  y::Array{Float64,1})

fac = y[4,end]/x[4];

    plot(x, y,
          xaxis=:log,
          yaxis=:log,
          xlabel="dt", ylabel="Error",
          color="blue",
          marker=:circle,
          legend=nothing)
    plot!(x, fac*x,
          xaxis=:log,
          yaxis=:log,
          color="red",
          linestyle=:dash)

end

# ==============================================================================
@time x, y = main()



##
function test(x::Array{Float64,1}, y::Array{Float64,1})

    xi = 1

    x1 = x[xi:end]
    y1 = y[xi:end]

    p1 = polyfit(log.(x1), log.(y1), 1)
    c = coeffs(p1)

    println(c)

    return x1, y1
end

x1, y1 = test(x,y)
