using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Statistics
using Quadmath

gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function largets_replica_runs(h5id, runs)
    # Only include one run per volume with the largest number of N_replicas
    data = [[read(h5id[r], "Nt"), read(h5id[r], "Ns"), read(h5id[r], "N_replicas")] for r in runs]
    maxr = similar(runs)
    for i in eachindex(data)
        matches = findall(x -> x[1:2] == data[i][1:2], data)
        j = findmax(x -> data[x][3], matches)[2]
        maxr[i] = runs[matches[j]]
    end
    return unique(maxr)
end
function cumulants(h5dset, run)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(h5dset, run)
    β = range(start = minimum(a), stop = maximum(a), length = 100)
    return cumulants(h5dset, run, β)
end
function cumulants(h5dset, run, β)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(h5dset, run)
    repeats = size(a, 2)

    CV = zeros((length(β), repeats))
    BC = zeros((length(β), repeats))

    for i in 1:repeats
        for j in eachindex(β)
            logZ = LLRParsing.log_partition_function(a[:, i], S, β[j], BigFloat)
            EN = LLRParsing.all_energy_moments(S, a[:, i], β[j], 4, Float128, BigFloat, logZ)
            u4 = EN[4] / (6.0 * V)^4
            u2 = EN[2] / (6.0 * V)^2
            u1 = EN[1] / (6.0 * V)^1
            CV[j, i] = 6V * (u2 - u1^2)
            BC[j, i] = 1 - u4 / u2 / u2 / 3
        end
    end
    return β, CV, BC
end
function beta_extremal(β, obs; f = findmax)
    repeats = size(obs, 2)
    βmax0 = zeros(repeats)
    for i in 1:repeats
        ind = f(obs[:, i])[2]
        βmax0[i] = β[ind]
    end
    βmax = mean(βmax0)
    Δβmax = std(βmax0) / sqrt(repeats)
    return βmax, Δβmax
end
function critical_beta_cumulants(h5dset, r; N, eps, min_iter = 5, max_iter = 20, w = 20)
    a = first(LLRParsing._set_up_histogram(h5dset, r))
    min_a, max_a = minimum(a), maximum(a)
    β = range(start = min_a, stop = max_a, length = N)

    βc_CV_old, βc_BC_old = +Inf, +Inf
    βc_CV, Δβc_CV = +Inf, +Inf
    βc_BC, Δβc_BC = +Inf, +Inf
    for i in 1:max_iter
        β, CV0, BC0 = cumulants(h5dset, r, β)

        βc_CV, Δβc_CV = beta_extremal(β, CV0; f = findmax)
        βc_BC, Δβc_BC = beta_extremal(β, BC0; f = findmin)
        βmin, βmax = extrema(β)
        βmin = min(βmin, βc_CV - w * Δβc_CV, βc_BC - w * Δβc_BC)
        βmax = max(βmin, βc_CV + w * Δβc_CV, βc_BC + w * Δβc_BC)
        βc_avg = (βc_CV + βc_BC) / 2
        β = range(start = (βmin + βc_avg) / 2, stop = (βmax + βc_avg) / 2, length = N)
        diff = max(abs(βc_CV_old - βc_CV), abs(βc_BC_old - βc_BC))
        βc_CV_old = βc_CV
        βc_BC_old = βc_BC
        if diff < eps && i > min_iter
            return βc_CV, Δβc_CV, βc_BC, Δβc_BC
        end
    end
    return βc_CV, Δβc_CV, βc_BC, Δβc_BC
end
function critical_beta_binder_cumulant(h5dset, r; N, eps, min_iter = 5, max_iter = 20, w = 20)
    a = first(LLRParsing._set_up_histogram(h5dset, r))
    min_a, max_a = minimum(a), maximum(a)
    β = range(start = min_a, stop = max_a, length = N)

    βc_BC_old = +Inf
    βc_BC, Δβc_BC = +Inf, +Inf
    for i in 1:max_iter
        β, CV0, BC0 = cumulants(h5dset, r, β)
        βc_BC, Δβc_BC = beta_extremal(β, BC0; f = findmin)
        βmin, βmax = extrema(β)
        βmin = min(βmin, βc_BC - w * Δβc_BC)
        βmax = max(βmin, βc_BC + w * Δβc_BC)
        β = range(start = (βmin + βc_BC) / 2, stop = (βmax + βc_BC) / 2, length = N)
        diff = abs(βc_BC_old - βc_BC)
        βc_BC_old = βc_BC
        if diff < eps && i > min_iter
            return βc_BC, Δβc_BC
        end
    end
    return βc_BC, Δβc_BC
end
function critical_beta_specific_heat(h5dset, r; N, eps, min_iter = 5, max_iter = 20, w = 20)
    a = first(LLRParsing._set_up_histogram(h5dset, r))
    min_a, max_a = minimum(a), maximum(a)
    β = range(start = min_a, stop = max_a, length = N)

    βc_CV_old = +Inf
    βc_CV, Δβc_CV = +Inf, +Inf
    for i in 1:max_iter
        β, CV0, BC0 = cumulants(h5dset, r, β)
        βc_CV, Δβc_CV = beta_extremal(β, CV0; f = findmax)
        βmin, βmax = extrema(β)
        βmin = min(βmin, βc_CV - w * Δβc_CV)
        βmax = max(βmin, βc_CV + w * Δβc_CV)
        β = range(start = (βmin + βc_CV) / 2, stop = (βmax + βc_CV) / 2, length = N)
        diff = abs(βc_CV_old - βc_CV)
        βc_CV_old = βc_CV
        if diff < eps && i > min_iter
            return βc_CV, Δβc_CV
        end
    end
    return βc_CV, Δβc_CV
end
function critical_cumulants_all_runs(h5file, outfile)
    h5dset = h5open(h5file)
    runs = keys(h5dset)
    runs = filter(!startswith("provenance"), runs)
    setprecision(BigFloat, 106)
    io = open(outfile, "w")
    print_provenance_csv(io)
    println(io, "group_family,Nc,n_replicas,run,T,L,βc_CV,Δβc_CV,βc_BC,Δβc_BC")
    for r in runs

        L = read(h5dset[r], "Ns")
        T = read(h5dset[r], "Nt")
        Nr = read(h5dset[r], "N_replicas")

        kws = (N = 50, eps = 1.0e-6)
        βc_CV, Δβc_CV = critical_beta_specific_heat(h5dset, r; kws...)
        βc_BC, Δβc_BC = critical_beta_binder_cumulant(h5dset, r; kws...)

        println(io, "Sp,4,$Nr,$r,$T,$L,$βc_CV,$Δβc_CV,$βc_BC,$Δβc_BC")
    end
    close(io)
    return nothing
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--outfile"
        help = "Where to save results"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    h5file = args["h5file"]
    outfile = args["outfile"]
    return critical_cumulants_all_runs(h5file, outfile)
end
main()
