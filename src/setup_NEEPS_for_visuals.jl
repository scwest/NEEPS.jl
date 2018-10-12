using NEEPS
days_to_event, event, clinical_patient_order = upload_clinical("../test/test1_clinical.txt")
expression_mat, element_order = upload_expression("../test/test1_expression.txt", clinical_patient_order)

null_size = 10000
min_threshold = 0.15
max_threshold = 0.85
num_workers = 4

addprocs(num_workers)
@everywhere using NEEPS

null_ps, lowest_pvals = parallel_null_and_curves(null_size, days_to_event,
event, min_threshold, max_threshold, expression_mat, num_workers)
