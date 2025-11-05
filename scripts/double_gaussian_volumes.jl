using LLRParsing
using HDF5
using Plots
using ArgParse
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
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
function plot_all_histogram_fits(file, plotfile, title)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    fid = h5open(file)
    runs = keys(fid)
    runs = filter(!startswith("provenance"), runs)
    runs = largets_replica_runs(fid, runs)
    plt = plot(title = title)
    for run in runs
        try
            βc = LLRParsing.beta_at_equal_heights(fid, run)
            plot_plaquette_histogram!(plt, fid, run, βc)
            plot!(plt, legend = :topright)
        catch
            @warn "Cannot determine critical β for $run"
        end
    end
    return savefig(plt, plotfile)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plotfile"
        help = "Where to save the plot"
        required = true
        "--title"
        help = "Title of the plot"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plotfile"]
    title = args["title"]
    return plot_all_histogram_fits(file, plotfile, title)
end
main()
