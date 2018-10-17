using Distributed
using SharedArrays

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

    null_jobs = SharedArray{Int, 1}((null_size))
    llp_jobs = SharedArray{Int, 1}((size(expression_mat)[1]))
    null_ps = SharedArray{Float64, 1}((null_size))
    lowest_pvals = SharedArray{Float64, 1}((size(expression_mat)[1]))
    directions = SharedArray{Float64, 1}((size(expression_mat)[1]))


    @everywhere include("survival_log_rank_pvals.jl")

    # create a channel with the jobs for the null distribution
    println("assigning parallel jobs")
    interval = 1
    null_left = length(null_jobs) - sum(null_jobs)
    llp_left = length(llp_jobs) - sum(llp_jobs)
    time_passed = 0
    io = IOBuffer()
    e = 0
    total = length(null_jobs) + length(llp_jobs)
    @sync begin
        #null_ps = zeros(null_size)
        for i in 1:null_size
            @spawn begin
                null_ps[i] = null_run(days_to_event, event, min_threshold, max_threshold)
                null_jobs[i] = 1
            end
            e += 1
            if e % 1000 == 0
                print("\r\tassigned $e of $total jobs")
                flush(io)
            end
        end

        #lowest_pvals = lowest_pvals = zeros(size(expression_mat)[1])
        for i in 1:size(expression_mat)[1]
            @spawn begin
                println("before llp")
                lowest_pvals[i], directions[i] = lowest_logrank_p(days_to_event, event, expression_mat[i,:], min_threshold, max_threshold)
                println("after llp")
                llp_jobs[i] = 1
                println("after llp jobs call")
            end
            e += 1
            if e % 1000 == 0
                print("\r\tassigned $e of $total jobs")
                flush(io)
            end
        end
        println("")

        println("running parallel jobs")
        println("\trefreshing every $interval seconds")
        while null_left > 0 || llp_left > 0
            null_left = length(null_jobs) - sum(null_jobs)
            llp_left = length(llp_jobs) - sum(llp_jobs)
            print("\r\tNull Jobs Left: $null_left\tSurvival Jobs Left: $llp_left\tTime Passed: $time_passed")
            flush(io)
            sleep(interval)
            time_passed += interval
        end
        println("")
    end

    println("all parallel jobs have finished")
    return null_ps, lowest_pvals, directions
end
