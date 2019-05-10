using PGFPlots

x = [1., 3., 4.5]
y = [1.5, 3., 1.]
r = 2
pts = Plots.Scatter(x, y, markSize=3);
c1 = Plots.Circle(x[1], y[1], r, style="black,very thick");
c2 = Plots.Circle(x[2], y[2], r, style="black,very thick");
c3 = Plots.Circle(x[3], y[3], r, style="black,very thick");
tri = Plots.Patch2D(hcat(x,y)', style="fill, fill opacity=0.3, draw opacity=0");
a = Axis([pts, c1, c2, c3, tri], xmin=-10,ymin=-10,xmax=10,ymax=10, hideAxis=true, width="7in", height="7in")
save("ceck.pdf", a)
