using Pkg
Pkg.activate(".")
using LLRParsing
using HDF5
using Plots
using MadrasSokal
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
    E0 = read(h5id, "$name/0/Rep_0/S0")

    E_therm = zeros(llr_therm, N_replicas, N_repeats)
    E_meas = zeros(llr_meas, N_replicas, N_repeats)
    E0 = zeros(N_replicas, N_repeats)

    for j in 0:(N_repeats - 1)
        for n in 0:(N_replicas - 1)
            E_meas[:, n + 1, j + 1] = read(h5id, "$name/$j/Rep_$n/E_meas")[:, end]
            E_therm[:, n + 1, j + 1] = read(h5id, "$name/$j/Rep_$n/E_therm")[:, end]
            E0[n + 1, j + 1] = read(h5id, "$name/$j/Rep_$n/S0")[end]
        end
    end

    # sort by central energy
    for r in axes(E_therm, 3)
        perm = sortperm(E0[:, r])
        E_therm[:, :, r] = E_therm[:, perm, r]
        E_meas[:, :, r] = E_meas[:, perm, r]
    end

    return E_therm, E_meas, ΔE, E0
end

function main(file, name;repeat = 1)
    h5id = h5open(file)
    E_therm, E_meas, ΔE, E0 = read_therm_meas(h5id, name)

    E_full = cat(E_therm,E_meas,dims=1)
    therms = size(E_therm,1) 
    
    plts = [ plot(E_full[:, i, repeat], ticks = :none, label = "", title = "replica #$i") for i in axes(E_meas, 2) ]
    plts = [ vspan!(plt, [1, therms], color = :blue, alpha = 0.2, labels = "therm") for plt in plts ]
    plt1 = plot(plts..., layout = grid(12, 8), size = (1000, 1500))
    
    plts = [ histogram(E_meas[:, i, repeat + 1], ticks = :none, label = "", title = "replica #$i") for i in axes(E_meas, 2) ]
    plt2 = plot(plts..., layout = grid(12, 8), size = (1000, 1500))
    plot!(plt2, plot_title = LLRParsing.fancy_title(name))

    savefig(plt1, "Trajectory.pdf")
    savefig(plt2, "Histogram.pdf")
    return nothing
end

file = "tmp/su4/su4.hdf5"
name = "5x32_96replicas"
main(file,name;repeat = 1)
