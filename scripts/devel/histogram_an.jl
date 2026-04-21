using Pkg
Pkg.activate(".")
using LLRParsing
using HDF5
using Plots
using LaTeXStrings
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

file = "data_assets/su4/all_su4_sorted.hdf5"
h5id = h5open(file)
runs = filter(startswith(r"[0-9]"), keys(h5id))
r = runs[5]
N_replicas = read(h5id[r], "N_replicas")

plts = []

for replica in 0:(N_replicas - 1)
    an = last.(LLRParsing.a_trajectory(h5id, r; replica))
    plt = histogram(an, label = L"a_n", ylabel = "", ticks = :none, title = "replica #$replica")
    push!(plts, plt)
end

pltF = plot(plts..., layout = grid(12, 8), size = (1000, 1500))
savefig(pltF, "Histogram_an.pdf")
