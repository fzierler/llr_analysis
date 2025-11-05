using ArgParse
using Plots
using LLRParsing
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
    left_margin = 1Plots.mm,
)

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
        "--critical_entropy"
        help = "CSV containing entropy at the critical point"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotfile = args["plot_file"]
    run = args["run_name"]
    entropy_csv = args["critical_entropy"]
    critical_entropy = Float64(readdlm(entropy_csv, skipstart = 1)[1, 1])
    return LLRParsing.plot_free_energy(file, plotfile, run, critical_entropy)
end
main()
