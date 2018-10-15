using NEEPS

# Obtain input from user and parse the clinical and expression files
clinical_patient_order, days_to_event, event, expression_mat,
element_order, min_threshold, max_threshold, null_size, output_prefix,
num_workers = get_input()

# Add workers if there are any to add
if num_workers > 1:
    addprocs(num_workers - nworkers())
    @everywhere using NEEPS
end

# Run parallel jobs:
# 1) Calculation of a null distribution using only the clinical samples
# 2) Calculation of the p-values for given expression patterns in exp. file
null_ps, lowest_pvals = parallel_null_and_curves(null_size, days_to_event,
event, min_threshold, max_threshold, expression_mat, num_workers)

# Construct output directory for both visuals and general
mkdir(output_prefix)
mkdir(string(output_prefix, "/visuals"))

# Generate null / lowest pval distributions early to clear up space
null_vs_lowest(null_ps, lowest_pvals,
string(output_prefix, "/visuals/null_vs_lowest.svg"))

# Run NEEP (This should be relatively quick)
ordered_neep_pvals = generate_neep_all(null_ps, lowest_pvals)

# Visuals for final NEEP values
alpha_choice(ordered_neep_pvals,
string(output_prefix, "/visuals/test1_alpha_choice.svg"))

# P-value Adjustment using MultipleTesting
ordered_neep_adj_pvals = adjust_neep_all(ordered_neep_pvals)

# Final p-value export
export_to_file(element_order, ordered_neep_adj_pvals,
string(output_prefix, "/neep_adjusted_pvalues.txt"))
