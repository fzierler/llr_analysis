using DelimitedFiles
using Plots
using ArgParse
using LaTeXStrings
gr(
    size = (350, 300),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 9,
    tickfontsize = 9,
    labelfontsize = 10,
    left_margin = 1Plots.mm,
)

function read_critical_betas(file, offset = 3)
    data, header = readdlm(file, ',', header = true, comments = true)
    Nt = data[:, 2 + offset]
    Ns = data[:, 3 + offset]
    A1 = data[:, 4 + offset]
    A2 = data[:, 5 + offset]
    βc = data[:, 6 + offset]
    Δβc = data[:, 7 + offset]
    return only(unique(Nt)), Ns, βc, Δβc, only(unique(A1)), only(unique(A2))
end
function plot_critical_beta!(plt, file; kws...)
    Nt, L, βc, Δβc, A1, A2 = read_critical_betas(file)
    tks = (inv.(L), (L"1/%$Li" for Li in L))
    lbl = L"$N_t=%$Nt$: ratio %$A1:%$A2"
    scatter!(plt, inv.(L), βc, xticks = tks, yerr = Δβc, label = lbl; kws...)
    return nothing
end
function plot_critical_beta(files, plotfile)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    plt = plot(ylabel = L"$\beta_{CV}(P)$", xlabel = L"$1/N_s$")
    markers = (:circle, :hexagon, :rect)
    colors = (:orange, :blue, :green)
    for (i, file) in enumerate(files)
        plot_critical_beta!(plt, file; markershape = markers[i], color = colors[i], markeralpha = 0.8)
    end
    plot!(plt; xflip = false, legend = :topright)
    return savefig(plt, plotfile)
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--plotfile"
        help = "Where to save the plot"
        required = true
        "arg"
        help = "CSV file(s) containing the values of critical beta"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plotfile = args["plotfile"]
    file = args["arg"]
    return plot_critical_beta(file, plotfile)
end
main()
