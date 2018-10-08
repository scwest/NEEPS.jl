function parallel_null_and_curves(null_size, days_to_event, event, min_threshold,
    max_threshold, expression_mat, num_workers)

    println("preparing parallel jobs")
    # moving existing variables to every process via forcing them not global
    let null_size = null_size
        @everywhere null_size = $null_size
    end
    let days_to_event = days_to_event
        @everywhere days_to_event = $days_to_event
    end
    let event = event
        @everywhere event = $event
    end
    let min_threshold = min_threshold
        @everywhere min_threshold = $min_threshold
    end
    let max_threshold = max_threshold
        @everywhere max_threshold = $max_threshold
    end
    let expression_mat = expression_mat
        @everywhere expression_mat = $expression_mat
    end

    @everywhere include("survival_log_rank_pvals.jl")
    include("transformation.jl")
    include("progress_bar.jl")

    # create a channel with the jobs for the null distribution
    println("running parallel jobs")
    @sync begin
        null_ps = zeros(null_size)
        for i in 1:null_size
            @spawn null_ps[i] = null_run(days_to_event, event, min_threshold, max_threshold)
        end

        lowest_pvals = lowest_pvals = zeros(size(expression_mat)[1])
        for i in 1:size(expression_mat)[1]
            @spawn lowest_pvals[i] = lowest_logrank_p(days_to_event, event, expression_mat[i], min_threshold, max_threshold)
        end
    end

    println("all parallel jobs have finished")
    return null_ps, lowest_pvals
end
