using ClusterComplex
using Plots
upscale = 1.0
default(size=(400*upscale,400*upscale))
pgfplots()
# gr()

X, H = ClusterComplex.dataset("TwoMoons")
s1 = scatter(X[1,:], X[2,:], 1, leg=false, markersize=1.7);
Plots.pdf(s1, "two-moons-ds")

X, H = ClusterComplex.dataset("Circles")
s2 = scatter(X[1,:], X[2,:], 1, leg=false, markersize=1.7);
Plots.pdf(s2, "circles-ds")

X, H = ClusterComplex.dataset("Sphere")
s3 = plot(X[1,:], X[2,:], X[3,:], seriestype=:scatter3d, legend=:none, markersize=1.7);
Plots.pdf(s3, "sphere-ds")

X, H = ClusterComplex.dataset("OIP300")
S = (rand(1:15000, 1200) |> unique)[1:1000]
s4 = plot(X[1,S], X[2,S], X[3,S], seriestype=:scatter3d, camera=(45, 15), legend=:none, markersize=1);
Plots.pdf(s4, "oip-300-ds")

X, H = ClusterComplex.dataset("OIP15")
S = (rand(1:15000, 1200) |> unique)[1:1000]
s5 = plot(X[1,S], X[2,S], X[3,S], seriestype=:scatter3d, camera=(45, 15), legend=:none, markersize=1);
Plots.pdf(s5, "oip-15-ds")
