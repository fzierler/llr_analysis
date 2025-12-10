using LLRParsing
using LaTeXStrings
using DelimitedFiles
using ArgParse
using Plots
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

function read_critical_betas(file, Nt0; offset = 0)
    data, header = readdlm(file, ',', header = true, comments = true)
    runs = data[:, 4]
    Nt = data[:, 5]
    Ns = data[:, 6]
    βc = data[:, 9 + offset]
    Δβc = data[:, 10 + offset]
    ind = findall(isequal(Nt0), Nt)
    return runs[ind], Nt[ind], Ns[ind], βc[ind], Δβc[ind]
end
function plot_critical_beta!(plt, Ns, βc, Δβc; kws...)
    tks = (inv.(Ns), (L"1/%$Li" for Li in Ns))
    scatter!(plt, inv.(Ns), βc, xticks = tks, yerr = Δβc; kws...)
    return nothing
end
function main_plot(file1, file2, outfile, Nt)
    runs11, Nt11, Ns11, βc11, Δβc11 = read_critical_betas(file1, Nt)
    runsCV, NtCV, NsCV, βcCV, ΔβcCV = read_critical_betas(file2, Nt; offset = -2)
    runsBC, NtBC, NsBC, βcBC, ΔβcBC = read_critical_betas(file2, Nt)

    @assert Nt11 == NtCV == NtBC
    Nt = first(Nt11)
    plt = plot(title = L"N_t = %$Nt", ylabel = L"$\beta_{CV}$", xlabel = L"$1/N_s$", legend = :topleft)

    plot_critical_beta!(plt, Ns11, βc11, Δβc11; ma = 0.7, marker = :circ, label = L"$\beta_{CV }(P)$")
    plot_critical_beta!(plt, NsBC, βcBC, ΔβcBC; ma = 0.7, marker = :rect, label = L"$\beta_{CV }(B_C)$")
    plot_critical_beta!(plt, NsCV, βcCV, ΔβcCV; ma = 0.7, marker = :diamond, label = L"$\beta_{CV }(C_V)$")
    savefig(outfile)

    return nothing
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--input_histogram"
        help = "CSV file with critical beta from histogram"
        required = true
        "--input_cumulants"
        help = "CSV file with critical beta from cumulants"
        required = true
        "--plot_file"
        help = "Where to save the table"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
    end
    return parse_args(s)
end
args = parse_commandline()
Nt = args["Nt"]
file1 = args["input_histogram"]
file2 = args["input_cumulants"]
file_plot = args["plot_file"]
main_plot(file1, file2, file_plot, Nt)
