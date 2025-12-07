using LLRParsing
using Plots
using HDF5
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
    left_margin = 7Plots.mm,
)

function overview(
        h5dset,
        run;
        repeat_id = 1,
        replica_id = read(h5dset[run], "N_replicas") ÷ 2,
    )
    Nreplicas = read(h5dset[run], "N_replicas")
    plt1 = full_trajectory_plot(h5dset, run, repeat_id, replica_id, lens = false)
    plt2 = full_trajectory_plot(h5dset, run, repeat_id, Nreplicas - 1, lens = false)
    plt3 = full_trajectory_plot(h5dset, run, repeat_id, 1, lens = false)
    plot!(plt2, legend = :topleft)
    plt = plot(plt3, plt1, plt2, layout = grid(1, 3), size = (1400, 1000))
    return plt
end
function overview_plot(file, run, plotfile)
    h5dset = h5open(file)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    Δa0 = a_vs_central_action(h5dset, run)[2]
    ind = findmax(Δa0)[2] - 1
    repeats = parse.(Int, read(h5dset[run], "repeats"))
    plt = overview(h5dset, run, repeat_id = first(repeats), replica_id = ind)
    return savefig(plt, plotfile)
end
function parse_commandline_per_run()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plots"
        required = true
        "--run_name"
        help = "Which dataset in the HDF5 file to plot"
        required = true
    end
    return parse_args(s)
end
function main_per_run()
    args = parse_commandline_per_run()
    return overview_plot(args["h5file"], args["run_name"], args["plot_file"])
end
main_per_run()
