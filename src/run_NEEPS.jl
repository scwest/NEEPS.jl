using NEEPS
clinical_patient_order, days_to_event, event, expression_mat,
element_order, min_threshold, max_threshold, null_size, output_filename,
num_workers = get_input()

addprocs(num_workers)
@everywhere using NEEPS

null_ps, lowest_pvals = parallel_null_and_curves(null_size, days_to_event,
event, min_threshold, max_threshold, expression_mat, num_workers)

null_vs_lowest(null_ps, lowest_pvals, "../test/test1_null_vs_lowest.svg")

ordered_neep_adj_pvals = generate_neep_all(null_ps, lowest_pvals)

export_to_file(element_order, ordered_neep_adj_pvals, output_filename)
