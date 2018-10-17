# Giving the User an immediate response upon tool use
println("setting up julia")
using Distributed
using NEEPS

# Obtain input from user and parse the clinical and expression files
clinical_patient_order, days_to_event, event, expression_mat,
element_order, min_threshold, max_threshold, null_size, output_prefix,
num_workers = get_input()

# Add workers if there are any to add
if num_workers > 1
    println("adding workers")
    addprocs(num_workers)
    @everywhere using NEEPS
end

# Run parallel jobs:
# 1) Calculation of a null distribution using only the clinical samples
# 2) Calculation of the p-values for given expression patterns in exp. file
null_ps, lowest_pvals, directions = parallel_null_and_curves(null_size, days_to_event,
event, min_threshold, max_threshold, expression_mat, num_workers)

# Construct output directory for both visuals and general
println("generating output directory")
mkdir(output_prefix)
mkdir(string(output_prefix, "/visuals"))

# Generate null / lowest pval distributions early to clear up space
println("visuals stage 1 of 2")
null_vs_lowest(null_ps, lowest_pvals,
string(output_prefix, "/visuals/null_vs_lowest.svg"))

# Run NEEP (This should be relatively quick)
ordered_neep_pvals = generate_neep_all(null_ps, lowest_pvals)

# Visuals for final NEEP values
println("visuals stage 2 of 2")
alpha_choice(sort(ordered_neep_pvals),
string(output_prefix, "/visuals/test1_alpha_choice.svg"))

# P-value Adjustment using MultipleTesting
ordered_neep_adj_pvals = adjust_neep_all(ordered_neep_pvals)

# Final p-value export
export_to_file(element_order, ordered_neep_adj_pvals, directions,
string(output_prefix, "/neep_adjusted_pvalues.txt"))
