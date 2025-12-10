using DelimitedFiles
using Plots
using ArgParse
using LaTeXStrings
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 1Plots.mm,
)

function read_critical_betas(file, Nt0, offset = 3)
    data, header = readdlm(file, ',', header = true, comments = true)
    Nt = data[:, 2 + offset]
    Ns = data[:, 3 + offset]
    A1 = data[:, 4 + offset]
    A2 = data[:, 5 + offset]
    βc = data[:, 6 + offset]
    Δβc = data[:, 7 + offset]
    ind = findall(isequal(Nt0), Nt)
    return only(unique(Nt[ind])), Ns[ind], βc[ind], Δβc[ind], only(unique(A1[ind])), only(unique(A2[ind]))
end
function plot_critical_beta!(plt, file, Nt0; kws...)
    Nt, L, βc, Δβc, A1, A2 = read_critical_betas(file, Nt0)
    tks = (inv.(L), (L"1/%$Li" for Li in L))
    lbl = L"$N_t=%$Nt$: ratio %$A1:%$A2"
    scatter!(plt, inv.(L), βc, xticks = tks, yerr = Δβc, label = lbl; kws...)
    return nothing
end
function plot_critical_beta(files, plotfile, Nt)
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    plt = plot(ylabel = L"$\beta_{CV}(P)$", xlabel = L"$1/N_s$")
    markers = (:circle, :hexagon, :rect)
    colors = (:orange, :blue, :green)
    for (i, file) in enumerate(files)
        plot_critical_beta!(plt, file, Nt; markershape = markers[i], color = colors[i], markeralpha = 0.8)
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
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
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
    Nt = args["Nt"]
    return plot_critical_beta(file, plotfile, Nt)
end
main()
