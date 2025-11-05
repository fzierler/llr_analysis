using HDF5
using LLRParsing
using Plots
using LaTeXStrings
using PDFmerger
using ProgressMeter
using ArgParse
gr(
    size = (350, 300),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 9,
    tickfontsize = 9,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function plot_repeat(fid, run, repeat, Nrep; kws...)
    title = LLRParsing.fancy_title(run) * ", repeat #$repeat"
    a0 = hcat([read(fid[run], "$repeat/Rep_$i/a_sorted") for i in 0:(Nrep - 1)]...)
    NRpRM = first(size(a0))
    plt = plot(
        a0;
        ms = 1,
        legend = :outerright,
        title,
        label = "",
        xlabel = "updates (NR + RM = $NRpRM)",
        ylabel = L"a_n",
    )
    return plot!(plt; kws...)
end
function plot_all_an_trajectories(file, outfile, run; kws...)
    fid = h5open(file)

    repeats = read(fid[run], "repeats")
    Nrep = read(fid[run], "N_replicas")

    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    isfile(outfile) && rm(outfile)
    tmp_plt_pdf = tempname() * ".pdf"
    for r in repeats
        plt = plot_repeat(fid, run, r, Nrep; kws...)
        savefig(plt, tmp_plt_pdf)
        append_pdf!(outfile, tmp_plt_pdf, cleanup = true)
    end
    return close(fid)

end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plots"
        required = true
        "--run_name"
        help = "LLR run in the HDF5 file to analyse"
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    return plot_all_an_trajectories(args["h5file"], args["plot_file"], args["run_name"])
end
main()
