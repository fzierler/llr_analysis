using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Statistics
using Peaks
using Quadmath
using DelimitedFiles
gr(
    size = (350, 300),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 9,
    tickfontsize = 9,
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
function cumulant_plots(h5file, Nt, critical_values)
    fid = h5open(h5file)
    pltCV = plot(legend = :topleft, xlabel = L"\beta", ylabel = L"C_V(\beta)", title = L"N_t = %$Nt")
    pltBC = plot(legend = :bottomright, xlabel = L"\beta", ylabel = L"B_V(\beta) - 2/3", title = L"N_t = %$Nt")
    runs = filter(!startswith("provenance"), keys(fid))
    runs = largets_replica_runs(fid, runs)

    # first determine a good plotting range
    β_min, β_max = +Inf, -Inf
    for r in runs
        a, _, _, _ = a_vs_central_action(fid, r)
        p_ind = only(findmaxima(a, 5).indices)
        m_ind = only(findminima(a, 5).indices)
        δ = m_ind - p_ind
        β_min, β_max = a[p_ind - δ], a[m_ind + 2δ]
    end

    for r in runs
        β = range(start = β_min, stop = β_max, length = 100)
        β, CV0, BC0 = cumulants(fid, r, β)
        repeats = size(CV0, 2)

        CV = dropdims(mean(CV0, dims = 2), dims = 2)
        ΔCV = dropdims(std(CV0, dims = 2), dims = 2) ./ sqrt.(repeats)
        BC = dropdims(mean(BC0, dims = 2), dims = 2)
        ΔBC = dropdims(std(BC0, dims = 2), dims = 2) ./ sqrt.(repeats)

        plot!(pltCV, β, CV, ribbon = ΔCV, label = LLRParsing.fancy_title(r), lw = 2)
        plot!(pltBC, β, BC .- 2 / 3, ribbon = ΔBC, label = LLRParsing.fancy_title(r), lw = 2)
    end
    plot!(pltCV, ylims = (0, maximum(ylims(pltCV))))
    plot!(pltBC, ylims = (minimum(ylims(pltBC)), 0))
    return pltCV, pltBC
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file_specific_heat"
        help = "Where to save the plot"
        required = true
        "--plot_file_binder_cumulant"
        help = "Where to save the plot"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        required = true
        arg_type = Int
        "--critical_values"
        help = "CSV file containing the critical values of beta"
        default = ""
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    pltCV, pltBC = cumulant_plots(args["h5file"], args["Nt"], args["critical_values"])
    savefig(pltCV, args["plot_file_specific_heat"])
    savefig(pltBC, args["plot_file_binder_cumulant"])
    return nothing
end
main()
