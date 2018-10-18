using DataFrames
using Gadfly
using Cairo

function null_vs_lowest(null_ps, lowest_pvals, outfilename)
    println("producing null vs lowest plot")
    null_ps = null_ps[1:length(lowest_pvals)]
    lpname = Array{String, 1}(undef, length(lowest_pvals))
    lpname[:] .= "Lowest p-value"
    lpboth = hcat(lowest_pvals, lpname)
    nuname = Array{String, 1}(undef, length(null_ps))
    nuname[:] .= "Null"
    nuboth = hcat(null_ps, nuname)
    ar = vcat(lpboth, nuboth)
    df = DataFrame(ar)
    rename!(df, :x1 => :vals)
    rename!(df, :x2 => :dist)

    p = plot(df, x="vals", color="dist", Geom.density,
    Guide.xlabel("p-value"), Guide.ylabel("density"),
    Guide.title("Distributions: Lowest p-value vs Expected Null"))
    draw(PNG(outfilename, 21cm, 15cm), p)
end

function make_blank_alpha_array(unp, alp)
    unp = unp[1:length(unp)/10]
    rhs = Array{Float64, 1}(undef, length(unp))
    for i in 1:length(unp)
        rhs[i] = alp*i/length(unp)
    end
    df = DataFrame(hcat(unp, rhs))
    return df
end

function alpha_choice(unordered_neep_pvals, outfilename)
    rhs = make_blank_alpha_array(unordered_neep_pvals, 0.05)
    plot1 = plot(rhs, x="x2", y="x1", Geom.point,
        intercept=[0], slope=[1], Geom.abline(color="white", style=:dash),
        Guide.xlabel("a*i/N"), Guide.ylabel("pval_i"),
        Guide.title("alpha = 0.05"))
    rhs = make_blank_alpha_array(unordered_neep_pvals, 0.1)
    plot2 = plot(rhs, x="x2", y="x1", Geom.point,
        intercept=[0], slope=[1], Geom.abline(color="white", style=:dash),
        Guide.xlabel("a*i/N"), Guide.ylabel("pval_i"),
        Guide.title("alpha = 0.1"))
    rhs = make_blank_alpha_array(unordered_neep_pvals, 0.2)
    plot3 = plot(rhs, x="x2", y="x1", Geom.point,
        intercept=[0], slope=[1], Geom.abline(color="white", style=:dash),
        Guide.xlabel("a*i/N"), Guide.ylabel("pval_i"),
        Guide.title("alpha = 0.2"))
    rhs = make_blank_alpha_array(unordered_neep_pvals, 0.3)
    plot4 = plot(rhs, x="x2", y="x1", Geom.point,
        intercept=[0], slope=[1], Geom.abline(color="white", style=:dash),
        Guide.xlabel("a*i/N"), Guide.ylabel("pval_i"),
        Guide.title("alpha = 0.3"))
    p = hstack(plot1, plot2, plot3, plot4)
    draw(PNG(outfilename, 60cm, 15cm), p)
end
