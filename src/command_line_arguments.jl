using ArgParse

"""
functions for parsing specific input arguments
"""
function upload_expression(input_filename, clinical_patient_order)
    fsize = open(input_filename) do infile
        fsize = 0
        for line in eachline(infile)
            fsize += 1
        end
        fsize
    end

    expression_mat, element_order =
    open(input_filename) do infile
        patient_order = split(strip(readline(infile)), ",")[2:end]
        indx = sortperm(patient_order, by=i->findfirst(clinical_patient_order.==i))
        element_order = String[]
        expression_mat = Matrix(undef, fsize, length(patient_order))
        i = 1
        for line in eachline(infile)
            line = split(strip(line), ",")
            push!(element_order, line[1])
            expressions = [parse(Float64, x) for x in line[2:end]]
            expressions = expressions[indx]
            expression_mat[i,:] = expressions'
            i += 1
        end
        (expression_mat, element_order)
    end
    return expression_mat, element_order
end

function upload_clinical(input_filename)
    days_to_event, event, clinical_patient_order =
    open(input_filename) do infile
        days_to_event = Int[]
        event = Int[]
        clinical_patient_order = String[]
        for line in eachline(infile)
            line = split(strip(line), ",")
            push!(clinical_patient_order, line[1])
            push!(days_to_event, parse(Float64, line[2]))
            push!(event, parse(Float64, line[3]))
        end
        combined = sort(collect(zip(days_to_event, event, clinical_patient_order)), by=x->x[1])
        days_to_event = getindex.(combined, 1)
        event = getindex.(combined, 2)
        clinical_patient_order = getindex.(combined, 3)
        (days_to_event, event, clinical_patient_order)
    end
    return clinical_patient_order, days_to_event, event
end


"""
defining variables from ARGS
"""
function get_input()
    println("reading command line")
    parsed_args = parse_commandline()

    println("uploading clinical file")
    clinical_patient_order, days_to_event, event =
    upload_clinical(parsed_args["clinical_filename"])

    println("uploading expression file")
    expression_mat, element_order =
    upload_expression(parsed_args["expression_filename"], clinical_patient_order)

    println("setting standard variables")
    min_threshold = parsed_args["min_thresh"]
    max_threshold = parsed_args["max_thresh"]
    null_size = parsed_args["null_size"]
    num_workers = parsed_args["num_workers"]

    output_prefix = parsed_args["output_prefix"]

    return clinical_patient_order, days_to_event, event, expression_mat,
    element_order, min_threshold, max_threshold, null_size, output_prefix,
    num_workers
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--clinical_filename", "-c"
            help = "file: [patient],[days_to_event],[event occured (0,1)]\\n"
        "--expression_filename", "-e"
            help = "file: table of expression values [genes/isoforms] by [patients]"
        "--min_thresh", "-i"
            help = "minimum checked threshold as proportion of samples"
            arg_type = Float64
            default = 0.15
        "--max_thresh", "-f"
            help = "maximum checked threshold as proportion of samples"
            arg_type = Float64
            default = 0.85
        "--output_prefix", "-o"
            help = "full path to where the desired output will be placed (directory)"
        "--null_size", "-n"
            help = "length of sampled null distribution (impacts precision)"
            arg_type = Int64
            default = 100000
        "--num_workers", "-w"
            help = "number of cores that you would like used"
            arg_type = Int64
            default = 4
    end

    return parse_args(s)
end
