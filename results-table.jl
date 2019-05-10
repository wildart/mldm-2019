import TOML

results = TOML.parsefile("results-final.toml")

dss = ["TwoMoons", "Sphere", "Circles", "OIP300", "OIP15"]
cols = ["vietorisrips", "witness", "mahalonobis_lmc", "mahalonobis_km"]
parts = Dict(1=>"\\% success", 2=>"median relative dominance", 4=>"median number of cells")

f = open("results-table.tex", "w")

for ds in dss
    data = results[ds]
    res = """
% dataset $ds
    \\multicolumn{1}{l}{\\textbf{$ds}} & & & & \\multicolumn{1}{l}{} \\\\ \\cline{1-5}
"""

    for part in sort(collect(keys(parts)))
        vals = fill("", size(cols)...)
        for (m,r) in data
            i = findfirst(isequal(m), cols)
            if i !== nothing
                vals[i] = string(part == 4 ? round(Int, r[part]) : round(r[part], digits=2))
            end
        end

        res *= """
% result 1
    \\multicolumn{1}{r|}{$(parts[part])} &
    \\multicolumn{1}{c}{$(vals[1])} &
    \\multicolumn{1}{c}{$(vals[2])} &
    \\multicolumn{1}{c}{$(vals[3])} &
    \\multicolumn{1}{c|}{$(vals[4])} \\\\ \\cline{1-5}
"""
    end
    println(f, res)
end

close(f)
