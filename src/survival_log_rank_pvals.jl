using Distributions


"""
Simple Functions
"""
# calculate the stepwise e value for generating log-rank statistic
get_e(m1, m2, n1, n2) = m1 - (n1 * (m1+m2) / (n1+n2))

# stepwise variance
function get_v(m1, m2, n1, n2)
    n = m1+m2 # wiki
    K = n1
    N = n1+n2
    return(n * (K/N) * ((N-K)/N) * ((N-n)/(N-1)))
end

# calculate variance from sum and sum of squares, and total
get_var(s::Float64, ss::Float64, total::Int64) = (ss - (s^2 / total)) / (total - 1)

# return ordered list of indices to change from group = 0 to group = 1
# also initializes group array
function get_ordered_indices(expression, min_threshold, max_threshold)
    group = Array{Bool, 1}(undef, length(expression))
    group[:] .= 0
    lowest_0_index = convert(Int16, floor(length(expression)*min_threshold))
    largest_threshold_index = convert(Int16, floor(length(expression)*max_threshold))
    p = sortperm(expression)[lowest_0_index:largest_threshold_index]
    group[p] .= 1
    p, group
end

"""
Survival Functions
"""
function lowest_logrank_p(days_to_event, event, expression,
    min_threshold, max_threshold)
    # days_to_event is sorted
    # indices of event correspond to indices of days_to_event
    # indices of expression correspond to indices of days_to_event
    # min/max threshold are for splitting KM curves
    # prepare for iterations
    p, group = get_ordered_indices(expression, min_threshold, max_threshold)
    # reverse p since pop is faster than shift
    p = reverse(p, 1)
    lowest_pval = 1.0
    direction = true
    while !isempty(p)
        #println(p)
        # move to the next threshold
        # find the index of the next lowest expression value
        group[pop!(p)] = 0  # set that spot to equal 0 (be in KM curve 'low')
        test_statistic, ndirection = get_test_statistic(days_to_event, event, group)
        pval = ccdf(Chisq(1), test_statistic)
        lowest_pval, direction = ifelse(pval <= lowest_pval, (pval, ndirection), (lowest_pval, direction))
    end
    return lowest_pval, direction
end

function get_test_statistic(days_to_event, event, group)
    n = [length(group)-sum(group), sum(group)]
    censored = [0, 0]
    observed = [0, 0]

    num = 0.0
    den = 0.0

    prev_days = days_to_event[1]
    for i in 1:length(days_to_event)
        if days_to_event[i] != prev_days
            num += get_e(observed[1], observed[2], n[1], n[2])
            den += get_v(observed[1], observed[2], n[1], n[2])

            n -= censored+observed
            censored = [0,0]
            observed = [0,0]
        end

        if event[i] == 0
            censored[group[i]+1] += 1
        else
            observed[group[i]+1] += 1
        end
        prev_days = days_to_event[i]
    end

    num += get_e(observed[1], observed[2], n[1], n[2])
    den += get_v(observed[1], observed[2], n[1], n[2])

    return (num^2)/den, num > 0
end

function null_run(days_to_event, event, min_threshold,
    max_threshold)
    expression = rand(length(days_to_event))
    llp, direction = lowest_logrank_p(days_to_event, event, expression, min_threshold,
    max_threshold)
    return llp
end
