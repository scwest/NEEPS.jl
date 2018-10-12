using MultipleTesting

"""
transformation of lowest-pvalues using estimated null distribution
"""
function generate_neep_pvals(sorted_lowest_pvals, null_ps)
    neep_pvals = zeros(length(sorted_lowest_pvals))
    neep_pvals[:] = 1
    pval = pop!(sorted_lowest_pvals) # lowest p-value pulled off
    spot = 1
    b = false
    for i in 1:length(null_ps)
        while pval < null_ps[i]
            neep_pvals[spot] = (i-1) / length(null_ps)
            spot += 1
            if !isempty(sorted_lowest_pvals)
                pval = pop!(sorted_lowest_pvals)
            else
                b = true
                break
            end
        end
        if b
            break
        end
    end
    return neep_pvals
end

function generate_neep_all(null_ps, lowest_pvals)
    # generate the NEEP p-values using the null distribution
    println("generating empirically estimated p-values")
    null_ps = sort(null_ps)
    println(lowest_pvals)
    sorted_lowest_pvals = sort(lowest_pvals, rev=true)
    indexed_lowest_pvals = sortperm(lowest_pvals)
    println(sorted_lowest_pvals)
    unordered_neep_pvals = generate_neep_pvals(sorted_lowest_pvals, null_ps) # from transformation.jl
    println(unordered_neep_pvals)
    ordered_neep_pvals = getindex.(sort(collect(zip(indexed_lowest_pvals, unordered_neep_pvals)), by=x->x[1]), 2)
    println(ordered_neep_pvals)
    ordered_neep_adj_pvals = adjust(ordered_neep_pvals, BenjaminiHochberg())
    return ordered_neep_adj_pvals
end
