using LLRParsing
using HDF5
using Statistics
using Plots
using LaTeXStrings
using LinearAlgebra
gr(
    size = (350, 300),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 9,
    tickfontsize = 9,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

function plot_pearson_correlation_matrix(ups, cor; title)
    ylabel = L"u_p"
    xlabel = L"u_p"
    plt = plot(; title, xlabel, ylabel, yflip = true)
    heatmap!(plt, ups, ups, abs.(cor), clims = (0, 1))
    return plt
end

file = "data_assets/Sp4_Nt5_sorted.hdf5"
fid = h5open(file)
typeof(fid)
runs = keys(fid)

for run in runs
    β0, βmin, βmax = 7.49, 7.488, 7.492
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(fid, run)
    beta = LLRParsing.beta_at_equal_heights(fid, run, β0, βmin, βmax)
    ups, prob, V, dS = LLRParsing.probability_density_repeats(a, S, beta, V)
    corP = cor(prob, dims = 2)
    covP = cov(prob, dims = 2) ./ size(prob)[2]
    cora = cor(a, dims = 2)
    P = dropdims(mean(prob, dims = 2), dims = 2)
    r = LLRParsing.fitting_range_double_gaussian(P; w = 5)

    title_a = "Pearson correlation matrix of \$a_n:\$" * LLRParsing.fancy_title(run)
    title_a = replace(title_a, "\$\$" => "")
    title_P = "Pearson correlation matrix of \$P_\\beta:\$" * LLRParsing.fancy_title(run)
    title_P = replace(title_P, "\$\$" => "")

    plt_a = plot_pearson_correlation_matrix(ups, cora; title = title_a)
    plt_P = plot_pearson_correlation_matrix(ups[r], corP[r, r]; title = title_P)
    #savefig(plt_a,"Pearson_correlation_an.pdf")
    #savefig(plt_P,"Pearson_correlation_P.pdf")
end
