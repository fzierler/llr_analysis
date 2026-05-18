using Pkg
Pkg.activate(".")
using LLRParsing
using HDF5
using Plots
using ArgParse
using Statistics
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

h5file = "data_assets/sp4/polyakov/polyakov_loop_data_4x20_64replicas.hdf5"
h5id = h5open(h5file)
ens = only(keys(h5id))

rx = r"(?<Nt>[0-9])x(?<Ns>[0-9]+)_(?<Nrep>[0-9]+)replicas"
m = match(rx, ens)
Nt, Ns, Nrep = parse(Int,m["Nt"]), parse(Int,m["Ns"]), parse(Int,m["Nrep"])
V = Nt*Ns^3

lp_abs  = read(h5id[ens],"lp_abs_samples")
lp_abs2 = read(h5id[ens],"lp_abs2_samples")

χlp_samples = @. (lp_abs2 - lp_abs^2)/V
χlp  = dropdims(mean(χlp_samples,dims=2),dims=2)
Δχlp = dropdims(std(χlp_samples,dims=2) ./ sqrt(size(χlp_samples,2)),dims=2)

plot(χlp, ribbon = Δχlp, label = LLRParsing.fancy_title(ens))