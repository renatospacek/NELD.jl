module Integrators

using Setup
using CellLists
using KR
using Compute

export EM!,
       BAO!,
       SEB!,
       SEAC!,
       SOILEB!, SOILEA!

# ==============================================================================
#                            Integrators
# ==============================================================================
function EM!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

    computeForce(X, Y, Z, λ)
    X.qtmp .= X.q .+ h.*X.p
    X.p += X.f.*h .+ λ.A*X.p.*h .- X.p.*h.*λ.γ .+ λ.A*X.q.*h.*λ.γ .+ X.dW3
    X.q .= X.qtmp
    X.prel .= X.p .- λ.A*X.q
    #addKE(X)
end

function BAO!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

    C1 = zeros(Float64, 3, 3)
    for i in 1:3
        C1[i,i] = exp(λ.A[i,i]*h)
    end

    X.p += 0.5*X.f*h
    X.q += h*X.p
    X.prel = X.p - λ.A*X.q

    PBC!(X, Y, λ)
    computeForce(X, Y, Z, λ)

    decay = exp(-λ.γ*h)
    diffuse = sqrt((1 - decay^2)/λ.β)

    X.p += 0.5*X.f*h
    X.p = C1*X.p

    X.p = X.p*decay + λ.A*X.q*(1 - decay) + X.dW
    X.prel = X.p - λ.A*X.q
end

function SEB!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

    computeForce(X, Y, Z, λ)
    for i in 1:3
        X.C1[i,i] = 1/(1 + λ.γ*h - λ.A[i,i]*h)
    end

    X.p = X.C1*(X.p + X.f*h + λ.A*X.q*h*λ.γ + X.dW3)
    X.q += h*X.p
    X.prel = X.p - λ.A*X.q
    addKE(X)
end

function SEAC!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

    X.q += h*X.p
    X.prel = X.p - λ.A*X.q

    Xtmp = Particle(nPart = λ.nPart)
    Xtmp.q .= X.q

    PBC!(Xtmp, Y, λ)
    computeForce(Xtmp, Y, Z, λ)
    X.f .= Xtmp.f

    addKE(X)

    X.p += λ.A*X.p*h - X.p*h*λ.γ + λ.A*X.q*h*λ.γ + X.dW3 + X.f*h
    X.prel = X.p - λ.A*X.q
end

function SOILEB!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

        X.p += (X.f + λ.A*X.q*λ.γ - λ.Γ2*X.p).*h.*0.5 - λ.Γ2*((X.f + λ.A*X.q*λ.γ - λ.Γ2*X.p).*h^2*0.5 + (X.dW4 - X.dW3*h*0.25).*2).*0.25 + X.dW3.*0.5;

     	X.q .= X.q .+ X.p.*h .+ X.dW4 - X.dW3.*h.*0.5;

     	X.prel .= X.p - λ.A*X.q;
     	Xtmp = Particle(nPart = λ.nPart);
        Xtmp.q .= X.q

        PBC!(Xtmp, Y, λ)
        computeForce(Xtmp, Y, Z, λ)
     	X.f .= Xtmp.f;

     	X.p += (X.f + λ.A*X.q*λ.γ - λ.Γ2*X.p).*h.*0.5 - λ.Γ2*((X.f + λ.A*X.q*λ.γ - λ.Γ2*X.p).*h^2*0.5  + (X.dW4 - X.dW3*h*0.25).*2).*0.25 + X.dW3.*0.5;

     	X.prel = X.p - λ.A*X.q;
end

function SOILEA!(X::Particle,
             Y::GenKR,
             Z::Clist,
             λ::Params,
             h::Float64)

        X.C2 .= 0.5.*(X.f + λ.A*X.q*λ.γ - λ.Γ*X.p).*h^2 .+ X.dW2

        X.Ftmp .= (X.f .+ λ.A*X.q.*λ.γ)

        X.q += X.p.*h .+ X.C2
        X.prel .= X.p .- λ.A*X.q;
        Xtmp = Particle(nPart = λ.nPart);
        Xtmp.q .= X.q

        PBC!(Xtmp, Y, λ)
        computeForce(Xtmp, Y, Z, λ)
        X.f .= Xtmp.f;

        #addKE(X)

        X.p += (X.Ftmp + X.f + λ.A*X.q.*λ.γ).*h.*0.5 - λ.Γ*(X.p.*h + X.C2) .+ X.dW1
        X.prel .= X.p .- λ.A*X.q
        addKE(X)
end

end
