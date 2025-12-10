using ArgParse
using Plots
using LLRParsing
using DelimitedFiles
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 1Plots.mm,
)

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plots"
        required = true
        "--critical_entropy"
        help = "CSV containing entropy at the critical point"
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
    Nt = args["Nt"]
    file = args["h5file"]
    plotfile = args["plot_file"]
    entropy_csv = args["critical_entropy"]
    critical_entropy = Float64(readdlm(entropy_csv, skipstart = 1)[1, 1])
    return LLRParsing.plot_entropy(file, plotfile, critical_entropy, Nt)
end
main()
