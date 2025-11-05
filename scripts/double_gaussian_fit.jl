using LLRParsing
using HDF5
using Plots
using ArgParse
using Peaks
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

function plot_all_histogram_fits(file, plotfile, run; fit, A1 = 1, A2 = 1, name = "fit")
    plotdir = dirname(plotfile)
    ispath(plotdir) || mkpath(plotdir)
    fid = h5open(file)
    plt = plot(title = LLRParsing.fancy_title(run))
    try
        βc = LLRParsing.beta_at_equal_heights(fid, run; A1, A2)
        plot_plaquette_histogram!(plt, fid, run, βc)
        if fit
            ups, f, Δf = LLRParsing.histogram_jackknife_fit(fid, run, βc)
            LLRParsing.plot_double_gaussian_fit_difference(plt, fid, run, βc)
            plot!(plt, ups, f, ribbon = Δf, label = "double Gaussian fit", lw = 2)
            # shade the region used in the fit
            ups, P, ΔP, covP, V, dS = probability_density(fid, run, βc)
            pks = findmaxima(P, 5)
            vspan!(plt, [ups[1], ups[pks.indices[1]]], color = :grey, alpha = 0.35, label = "")
            vspan!(plt, [ups[pks.indices[2]], ups[end]], color = :grey, alpha = 0.35, label = "")
        end
    catch
        @warn "Error for plot with peak-ratio $A1:$A2 for $run"
    end
    savefig(plotfile)
    return close(fid)
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
        "--run_name"
        help = "Run to be analysed"
        required = true
        "--peak1"
        help = "Ratio of peak heights used for finding the critical coupling"
        arg_type = Int
        default = 1
        "--peak2"
        help = "Ratio of peak heights used for finding the critical coupling"
        arg_type = Int
        default = 1
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plot_file"]
    run_name = args["run_name"]
    A1 = args["peak1"]
    A2 = args["peak2"]
    plot_all_histogram_fits(file, plotfile, run_name; fit = true, A1, A2, name = "")
    return nothing
end
main()
