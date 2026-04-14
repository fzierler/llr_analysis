using LLRParsing
using HDF5
using Plots
gr(
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 5,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
    top_margin = 0Plots.mm,
)
function read_therm_meas(h5id, name)
    N_repeats = read(h5id, "$name/N_repeats")
    N_replicas = read(h5id, "$name/N_replicas")
    llr_therm = read(h5id, "$name/0/Rep_0/llr_therm")
    llr_meas = read(h5id, "$name/0/Rep_0/llr_meas")
    ΔE = read(h5id, "$name/0/Rep_0/dS0")

    E_therm = zeros(llr_therm, N_replicas, N_repeats)
    E_meas = zeros(llr_meas, N_replicas, N_repeats)

    for j in 0:(N_repeats - 1)
        for n in 0:(N_replicas - 1)
            E_meas[:, n + 1, j + 1] = read(h5id, "$name/$j/Rep_$n/E_meas")[:, end]
            E_therm[:, n + 1, j + 1] = read(h5id, "$name/$j/Rep_$n/E_therm")[:, end]
        end
    end

    return E_therm, E_meas, ΔE
end

h5id = h5open("tmp/su4/su4.hdf5")
name = "5x32_96replicas"

N_repeats = read(h5id, "$name/N_repeats")
N_replicas = read(h5id, "$name/N_replicas")
Nt = read(h5id, "$name/Nt")
Ns = read(h5id, "$name/Ns")

E_therm, E_meas, ΔE = read_therm_meas(h5id, name)

repeat = 0
plts = [ histogram(E_meas[:, i, repeat + 1], ticks = :none, label = "", title = "replica #$i") for i in axes(E_meas, 2) ]
plt = plot(plts..., layout = grid(12, 8), size = (1000, 1500))
plot!(plt, plot_title = LLRParsing.fancy_title(name))
savefig("histograms.pdf")
