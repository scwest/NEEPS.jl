using Distributions


"""
Simple Functions
"""
# calculate the stepwise e value for generating log-rank statistic
get_e(m1, m2, n1, n2) = n1 * (m1+m2) / (n1+n2)

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
    println("before isempty")
    while !isempty(p)
        #println(p)
        # move to the next threshold
        # find the index of the next lowest expression value
        group[pop!(p)] = 0  # set that spot to equal 0 (be in KM curve 'low')
        println("before gtt")
        test_statistic, ndirection = get_test_statistic(days_to_event, event, group)
        println("after gtt")
        pval = ccdf(Chisq(1), test_statistic)
        println("cant be that")
        lowest_pval, direction = ifelse(pval <= lowest_pval, (pval, ndirection), (lowest_pval, direction))
        println("ugh")
    end
    return lowest_pval, direction
end

function get_test_statistic(days_to_event, event, group)
    total = [0,0]
    n1s = sum(group)
    n = [length(group)-n1s, n1s]
    m = [0,0]
    g = [0,0]
    s = 0.0
    ss = 0.0
    prev_days = 0
    for i in 1:length(days_to_event)
        if days_to_event[i] == prev_days || sum(m) == 0
            if event[i] == 0
                g[group[i]+1] += 1
            else
                n -= g
                g = [0,0]
                m[group[i]+1] += 1
            end
        else
            e = get_e(m[1], m[2], n[1], n[2])
            s += e
            prev_days = days_to_event[i]

            total += m
            n -= m+g
            m = [0,0]
            g = [0,0]
            if event[i] == 0
                g[group[i]+1] += 1
            else
                m[group[i]+1] += 1
            end
        end
    end
    a = (total[1] - s)^2 / s
    s2 = total[1]+total[2]-s
    b = (total[2] - s2)^2 / s2

    return a+b, (total[1] - s) > 0
end

function null_run(days_to_event, event, min_threshold,
    max_threshold)
    expression = rand(length(days_to_event))
    llp, direction = lowest_logrank_p(days_to_event, event, expression, min_threshold,
    max_threshold)
    return llp
end
