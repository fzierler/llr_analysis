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

function plot_all_histogram_fits(file, plotfile, title, Nt)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    fid = h5open(file)
    runs = filter(!startswith("provenance"), keys(fid))
    runs = filter(r -> read(fid[r], "Nt") == Nt, runs)
    plt = plot(title = title)
    for run in runs
        try
            βc = LLRParsing.beta_at_equal_heights(fid, run)
            plot_plaquette_histogram!(plt, fid, run, βc)
            plot!(plt, legend = :outerright)
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
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        required = true
        arg_type = Int
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plotfile"]
    title = args["title"]
    Nt = args["Nt"]
    return plot_all_histogram_fits(file, plotfile, title, Nt)
end
main()
