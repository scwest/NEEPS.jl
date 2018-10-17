function export_to_file(element_order, ordered_neep_adj_pvals, directions, output_filename)
    print("exporting p-values")
    open(output_filename, "w") do outfile
        for (element, pval, direction) in zip(element_order, ordered_neep_adj_pvals, directions)
            dstring = ifelse(direction > 0, "high expression survived longer", "low expression survived longer")
            write(outfile, "$element\t$pval\t$dstring\n")
        end
    end
end
