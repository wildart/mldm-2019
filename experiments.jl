using ClusterComplex
import ComputationalHomology: ComputationalHomology, filtration
import Clustering
import LMCLUS
import TOML
import Statistics: median, mean
import Random

const SEED = 921369837
Random.seed!(SEED)
prngs = [Random.MersenneTwister(SEED)]

KMeansPart = (Y, k)->Clustering.kmeans(collect(Y), k)
LMCLUSPart = (Y, bb)->let p = LMCLUS.Parameters(size(Y, 1))
           p.sampling_heuristic=2
           p.best_bound = abs(bb)
           p.basis_alignment = true
           p.bounded_cluster = true
           p.min_cluster_size = 30
           p.sampling_factor = 0.1
           res = LMCLUS.LMCLUSResult(LMCLUS.lmclus(Y, p, prngs)...)
           if bb < 0
                dfun = (X,m)  -> LMCLUS.distance_to_manifold(X, mean(m), LMCLUS.projection(m))
                efun = (X,ms) -> LMCLUS.MDL.calculate(LMCLUS.MDL.OptimalQuant, ms, X, 32, 16)
                res = LMCLUS.refine(res, Y, dfun, efun, bounds=true)
           end
           res
end

constractions = Dict(
    "vietorisrips" => (X, r)->ComputationalHomology.vietorisrips(X, r, maxoutdim=size(X,1)-1),
    "witness" => (X, r, l, hom=size(X,1)-1)->ComputationalHomology.witness(X, l, r, maxoutdim=hom),
    "lmcdist_km" => (X, r, k, hom=size(X,1)-1)->clustercomplex(X, KMeansPart(X, k), r, maxoutdim=hom, expansion = :inductive, method=:subspacemahalonobis),
    "lmcdist_lmc" => (X, r, bb, hom=size(X,1)-1)->clustercomplex(X, LMCLUSPart(X, bb), r, maxoutdim=hom, expansion = :inductive, method=:subspacemahalonobis),
    "mahalonobis_km" => (X, r, k, hom=size(X,1)-1)->clustercomplex(X, KMeansPart(X, k), r, maxoutdim=hom, expansion = :inductive),
    "mahalonobis_lmc" => (X, r, bb, hom=size(X,1)-1)->clustercomplex(X, LMCLUSPart(X, bb), r, maxoutdim=hom, expansion = :inductive),
    "mahalonobis_km_nu1" => (X, r, k, hom=size(X,1)-1)->clustercomplex(X, KMeansPart(X, k), r, maxoutdim=hom, expansion = :inductive, ν=1),
    "mahalonobis_lmc_nu1" => (X, r, bb, hom=size(X,1)-1)->clustercomplex(X, LMCLUSPart(X, bb), r, maxoutdim=hom, expansion = :inductive,ν=1),
)

function prepareparams(params, pnames = ["r", "l", "bb", "k", "hom"])
    res = Any[]
    for k in pnames
        haskey(params, k) && push!(res, params[k])
    end
    return tuple(res...)
end

function experiment(datasets, mthds; N=100)
    constraction_params = TOML.parsefile("experiments.toml")
    results = Dict(ds => Dict(m => Any[] for m in mthds) for ds in datasets)

    for ds in datasets
        X, H, L = ClusterComplex.dataset(ds)
        println("$ds: $H")

        for m in mthds
            print("Run $m:")
            n = N
            if m == "vietorisrips"
                n = 1
                if ds == "OIP300" || ds == "OIP15"
                    n = 0
                end
            end
            rdom = zeros(n)
            clen = zeros(Int, n)
            psrt = zeros(n)
            for i in 1:n
                divs = get(constraction_params[ds][m], "divs", Inf)
                try
                    res = constractions[m](X, prepareparams(constraction_params[ds][m])...)
                    flt = filtration(res[1:2]...; divisions=divs)
                    RD, R0, _ = ClusterComplex.dominance(flt, H)
                    rdom[i] = RD
                    psrt[i] = R0
                    if !isnan(RD)
                        cplx = complex(flt, R0)
                        clen[i] = length(cplx)
                    end
                    print(isnan(RD) ? "x" : ".")
                catch err
                    print("$i")
                end
            end
            succ = n == 1 ? [1] : findall(rdom .> 0)
            push!(results[ds][m], 100*length(succ)/N)
            push!(results[ds][m], length(succ) > 0 ? median(rdom[succ]) : 0.0)
            push!(results[ds][m], length(succ) > 0 ? mean(rdom[succ]) : 0.0)
            push!(results[ds][m], length(succ) > 0 ? median(clen[succ]) : 0.0)
            push!(results[ds][m], length(succ) > 0 ? mean(clen[succ]) : 0.0)
            push!(results[ds][m], length(succ) > 0 ? mean(psrt[succ]) : 0.0)
            println("\nmedian relative dominance: $(results[ds][m][2]) [$(results[ds][m][1])%]")
        end
    end

    return results
end

datasets = ["Sphere","TwoMoons","Circles","OIP300", "OIP15"]
# datasets = ["OIP300", "OIP15"]
# datasets = ["OIP15"]

mthds = ["witness","lmcdist_lmc","mahalonobis_lmc","lmcdist_km","mahalonobis_km","vietorisrips"]
# mthds = ["lmcdist_km","mahalonobis_km"]
# mthds = ["lmcdist_lmc","mahalonobis_lmc"]
# mthds = ["mahalonobis_lmc","mahalonobis_km"]
# mthds = ["witness", "vietorisrips"]
# mthds = ["witness","mahalonobis_lmc","mahalonobis_km","vietorisrips"]
# mthds = ["witness"]
# mthds = ["vietorisrips"]

res = experiment(datasets, mthds)
open(io->TOML.print(io, res), "results.toml", "w")
mv("results.toml", "results-$(time_ns()).toml")
# success rate, median dominance, mean dominance, median complex size, mean complex size, profile start F-value
