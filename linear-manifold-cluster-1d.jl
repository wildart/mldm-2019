using LinearManifoldCluster
import Random
using Statistics
using MultivariateStats
using LinearAlgebra

import LaTeXStrings: latexstring
using Plots
upscale = 2
Plots.scalefontsizes(upscale)
default(size=(600*upscale,400*upscale))
pgfplots()
# gr()
# clibrary(:colorcet)
# theme(:wong2)

const SEED = 921369837
Random.seed!(SEED)

R(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

function elliptical_bound(S::AbstractMatrix{<:Real}, χ::Real)
    F = eigen(Symmetric(S))
    λ = F.values
    FI = sortperm(λ, rev=true)
    EV = F.vectors[:,FI[1]]
    ϕ = atan(EV[2], EV[1])
    if ϕ < 0
        ϕ += 2π
    end

    θ = range(0.0, 2π, step=0.1)
    a, b = χ*sqrt(λ[FI[1]]), χ*sqrt(λ[FI[2]])
    x, y  = a.*cos.(θ), b.*sin.(θ)
    ellipse = R(ϕ) * hcat(x, y)'
    return ellipse
end

# generate
n = 200          # number of points in generated manifold
N = 2            # space dimensionality
M = [1]         # number of manifolds
τ = .2 => .8   # translation vector bounds
E = [[.12]]      # bound on an extend from manifold for all manifolds)
Φ = Vector{Float64}[[.9]]   # bounds on manifold points (per manifold wrt dims)
X = generate(n,N,M,τ,Φ,E)[1][1]

# rotate
RX = R(π/5)*X
RX .+=  [0.5; -0.13]

# perform PCA
M = fit(PCA, RX)
P = projection(M)
O = mean(M)

e = 0.1
L1 = [O+P[:,1] O-P[:,1]]
L2 = [O+P[:,2]*e O-P[:,2]*e]

plot(RX[1,:], RX[2,:], seriestype=:scatter, aspect_ratio=:equal,
     xlims=(0.0,1.0), ylims=(0.2,0.85), markersize=2,
     legend=:bottomright, lab=latexstring("LM cluster, \$C_{\\Lambda,\\theta}\$"), palette=:viridis)
plot!(L1[1,:], L1[2,:], seriestype=:line, linewidth=2, lab=latexstring("linear manifold, \$\\Lambda\$"))
plot!(L2[1,:], L2[2,:], seriestype=:line, linewidth=2, linestyle=:dash, lab="orthoganal subspace")


# boundary
C = cov(RX.-O, dims=2)

for χ in [2.0, 2.5, 3.0]
    B = elliptical_bound(C, χ).+O
    plot!(B[1,:], B[2,:], seriestype=:path, linewidth=1,
          lab=latexstring("geom. LM cluster, \$C_{\\Lambda,\\chi=$χ}\$"))
end

Plots.pdf("linear-manifold-cluster-1d")
