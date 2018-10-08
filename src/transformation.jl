using MultipleTesting

"""
transformation of lowest-pvalues using estimated null distribution
"""
function generate_neep_pvals(sorted_lowest_pvals, null_ps)
    neep_pvals = zeros(length(sorted_lowest_pvals))
    pval = pop!(sorted_lowest_pvals)
    spot = 1
    for i in 1:length(null_ps)
        while pval < null_ps[i]
            neep_pvals[spot] = i / length(null_ps)
            spot += 1
            if !isempty(sorted_lowest_pvals)
                pval = pop!(sorted_lowest_pvals)
            else
                break
            end
        end
        if isempty(sorted_lowest_pvals)
            break
        end
    end

    if !isempty(sorted_lowest_pvals)
        for _ in 1:length(sorted_lowest_pvals)
            neep_pvals[spot] = 1
            spot += 1
        end
    end

    return neep_pvals
end

function generate_neep_all(null_ps, lowest_pvals)
    # generate the NEEP p-values using the null distribution
    print("generating empirically estimated p-values")
    null_ps = sorted(null_ps)
    sorted_lowest_pvals = sort(lowest_pvals, rev=true)
    indexed_lowest_pvals = sortperm(lowest_pvals, rev=true)
    unordered_neep_pvals = generate_neep_pvals(sorted_lowest_pvals, null_ps) # from transformation.jl
    ordered_neep_pvals = getindex.(sort(collect(zip(indexed_lowest_pvals, unordered_neep_pvals)), by=x->x[1]), 1)
    ordered_neep_adj_pvals = adjust(ordered_neep_pvals, BenjaminiHochberg())
    return ordered_neep_adj_pvals
end
