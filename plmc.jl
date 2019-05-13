using Plots
using TDA
using ClusterComplex
import ComputationalHomology: ComputationalHomology, filtration, intervals
import Random
import StatsBase

const SEED = 921369837

import Clustering
function kmeans(Y, k)
    Random.seed!(SEED)
    Clustering.kmeans(Y, k)
end

probmah(p) = sqrt(-2*log(1-p))

X, H, L = ClusterComplex.dataset("TwoMoons")

CL = kmeans(X, 7)
plot(CL, X, markersize=1.5, legend=:none, color_palette=:viridis)

α, αmax = 0.05, 0.001
res = clustercomplex(X, CL, probmah(1-αmax), maxoutdim=1)
plot!(res[3], probmah(1-α), colors=true, color_palette=:viridis)

flt = filtration(res[1:2]...)
plot!(complex(flt, probmah(1-α)), CL.centers',
    linewidth=3, color_palette=:viridis)

Plots.pdf("plmc")
