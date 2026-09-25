using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Peaks
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

a_vs_central_action_plot(h5id, runs::Vector; kws...) =
    a_vs_central_action_plot!(plot(), h5id, runs; kws...)
function largets_replica_runs(h5id, runs)
    # Only include one run per volume with the largest number of N_replicas
    data = [
        [read(h5id[r], "Nt"), read(h5id[r], "Ns"), read(h5id[r], "N_replicas")] for
            r in runs
    ]
    maxr = similar(runs)
    for i in eachindex(data)
        matches = findall(x -> x[1:2] == data[i][1:2], data)
        j = findmax(x -> data[x][3], matches)[2]
        maxr[i] = runs[matches[j]]
    end
    return unique(maxr)
end
function a_vs_central_action_plot!(plt, h5id, runs::Vector; kws...)
    # default plot limits
    xmin, xmax = +Inf, -Inf
    ymin, ymax = +Inf, -Inf
    for run in runs
        a0, Δa0, S0, _ = a_vs_central_action(h5id, run)
        Nt = read(h5id[run], "Nt")
        Ns = read(h5id[run], "Ns")
        Nrep = read(h5id[run], "N_replicas")
        up = @. S0 / (6 * Ns^3 * Nt)
        label = L"N_s\!=\!%$Ns,N_{\!\mathrm{rep}}\!=\!%$Nrep"
        LLRParsing.a_vs_central_action_plot!(plt, a0, Δa0, S0, Nt, Ns, Nrep; label, kws...)
        # find useful plot limits for the volume comparison
        p_ind = findmaxima(a0, 5).indices
        m_ind = findminima(a0, 5).indices
        if length(m_ind) == length(p_ind) == 1
            p_ind = only(p_ind)
            m_ind = only(m_ind)
            # y-axis range: Show peak/valley and some extra
            s = 0.25*(a0[p_ind] - a0[m_ind])
            ymin, ymax = min(ymin, a0[m_ind] - s ), max(ymax, a0[p_ind] + s )
            # x-axis range: Stop plot where a_n 'leaves' the y-axis range 
            x_min_ind = findfirst(i -> a0[i] > ymin , eachindex(a0))
            x_max_ind = findlast(i -> a0[i] < ymax , eachindex(a0))
            xmin, xmax = min(xmin, up[x_min_ind]), max(xmax, up[x_max_ind])
            plot!(plt, xlims = (xmin, xmax), ylims = (ymin, ymax))
        end
    end
    return plt
end
function an_action_volumes(file, plotdest, Nt, Ns, Nrep; title, largets_replicas, zoom)
    ispath(dirname(plotdest)) || mkpath(dirname(plotdest))
    h5id = h5open(file)
    runs = keys(h5id)
    runs = filter(!startswith("provenance"), runs)
    if !iszero(Nt)
        runs = filter(r -> read(h5id[r], "Nt") == Nt, runs)
    end
    if !iszero(Ns)
        runs = filter(r -> read(h5id[r], "Ns") == Ns, runs)
    end
    if !iszero(Nrep)
        runs = filter(r -> read(h5id[r], "N_replicas") == Nrep, runs)
    end
    if largets_replicas
        runs = largets_replica_runs(h5id, runs)
    end
    plt = a_vs_central_action_plot(h5id, runs, lens = false)
    if !isempty(title)
        title = latexstring(title)
    end
    plot!(plt; legend = :bottomright, xlabel = L"u_p", ylabel = L"a_n", title)
    return savefig(plt, plotdest)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plot"
        required = true
        "--title"
        help = "Title of the plot"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
        "--Ns"
        help = "Ns of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
        "--Nrep"
        help = "Nrep of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
        "--largets_replicas"
        help = "include only run with largest number of replicas"
        default = true
        arg_type = Bool
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotdst = args["plot_file"]
    title = args["title"]
    Nt = args["Nt"]
    Ns = args["Ns"]
    Nrep = args["Nrep"]
    largets_replicas = args["largets_replicas"]
    return an_action_volumes(file, plotdst, Nt, Ns, Nrep; title, largets_replicas)
end
main()
