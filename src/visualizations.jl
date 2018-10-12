using DataFrames
using StatPlots
using Plots
using GR

function null_vs_lowest(null_ps, lowest_pvals, outfilename)
    null_ps = null_ps[1:length(lowest_pvals)]
    lpname = Array{String, 1}(length(lowest_pvals))
    lpname[:] = "Lowest p-value"
    lpboth = hcat(lowest_pvals, lpname)
    nuname = Array{String, 1}(length(null_ps))
    nuname[:] = "Null"
    nuboth = hcat(null_ps, nuname)
    ar = vcat(lpboth, nuboth)
    df = DataFrame(ar)
    rename!(df, :x1 => :vals)
    rename!(df, :x2 => :dist)

    gr()
    beginprint(outfilename)
    @df df density(:vals, group=(:dist))
    plot!(xlabel="p-value", ylabel="density")
    plot!(fmt=:svg)
    endprint()

end
