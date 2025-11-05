using LLRParsing
using HDF5
using Plots
using LaTeXStrings
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

file = "data_assets/Sp4_Nt5_sorted.hdf5"
h5id = h5open(file)
runs = filter(startswith(r"[0-9]"), keys(h5id))
r = runs[8]
N_replicas = read(h5id[r], "N_replicas")

for replica in 1:N_replicas
    an = last.(LLRParsing.a_trajectory(h5id, r; replica))
    title = LLRParsing.fancy_title(r) * " - replica: $replica"
    plt = plot(; title, xlabel = L"a_n", ylabel = "counts")
    histogram!(plt, an, label = "")
    display(plt)
end
