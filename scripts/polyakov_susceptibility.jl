using Pkg
Pkg.activate(".")
using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Statistics
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 9,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function plot_polyakov_loop(h5file,plot_name,Nt_plot)
    plt = plot(title="Polyakov loop susceptibility",xlabel=L"\beta",ylabel=L"\chi_{|l_p|}/V")

    h5id = h5open(h5file)
    runs = keys(h5id)
    runs = filter(!startswith("provenance"), runs)

    for ens in runs
        
        rx = r"(?<Nt>[0-9])x(?<Ns>[0-9]+)_(?<Nrep>[0-9]+)replicas"
        m = match(rx, ens)
        Nt, Ns, Nrep = parse(Int,m["Nt"]), parse(Int,m["Ns"]), parse(Int,m["Nrep"])
        V = Nt*Ns^3

        if Nt == Nt_plot
            β = read(h5id[ens],"beta")
            lp_abs  = read(h5id[ens],"lp_abs_samples")
            lp_abs2 = read(h5id[ens],"lp_abs2_samples")

            χlp_samples = @. (lp_abs2 - lp_abs^2)/V
            χlp  = dropdims(mean(χlp_samples,dims=2),dims=2)
            Δχlp = dropdims(std(χlp_samples,dims=2) ./ sqrt(size(χlp_samples,2)),dims=2)

            plot!(plt, β, χlp, ribbon = Δχlp, label = LLRParsing.fancy_title(ens))
        end
    end

    savefig(plt,plot_name)
    return nothing
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--plotfile"
        help = "Where to write the resulting HDF5 file containing the sorted results"
        required = true
        "arg"
        help = "HDF% files to be included in the plot"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plotfile = args["plotfile"]
    h5files = args["arg"]
    plot_polyakov_loop(h5files,plotfile)
end
main()