function export_all_pvals(all_pvals_mat, all_thresholds, element_order, output_filename)
    open(output_filename, "w") do outfile
        write(outfile, string(",", join(all_thresholds, ","), "\n"))
        for i in 1:length(element_order)
            write(outfile, string(element_order[i], ",", join(all_pvals_mat[i,:], ","), "\n"))
        end
    end
end

function export_to_file(element_order, lowest_pvals, ordered_neep_pvals,
                        ordered_neep_adj_pvals, directions, threshs,
                        output_filename)
    println("exporting p-values")
    open(output_filename, "w") do outfile
        write(outfile, "name\tlowest_p-value\tneep_p-value\tadjusted_neep\tdirection\n")
        for (element, lp, neep, adj, direction, thresh) in
            zip(element_order, lowest_pvals, ordered_neep_pvals, ordered_neep_adj_pvals, directions, threshs)
            dstring = ifelse(direction > 0, "high expression survived longer", "low expression survived longer")
            write(outfile, "$element\t$lp\t$neep\t$adj\t$dstring\t$thresh\n")
        end
    end
end

function export_distribution(nps, output_filename)
    open(output_filename, "w") do outfile
        for npval in nps
            write(outfile, "$npval\n")
        end
    end
end
