using Pkg
Pkg.activate(".")
using LaTeXStrings
using LLRParsing
using Clustering
using Plots
using HDF5
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

h5dset = h5open("data_assets/su4/all_su4_sorted.hdf5")
run = "5x32_112replicas_run2"

h5dset = h5open("data_assets/sp4/all_sp4_sorted.hdf5")
run = "4x28_64replicas_run2"

# read last values of an
a, S, Nt, Ns, V = LLRParsing._set_up_histogram(h5dset, run)
up = S/(6V)

# pre-allocate an array for the results 
means = zeros(size(a,1),2)
for i in axes(a,1) 
    # change layout to the from expected by Clustering.jl 
    # (it requires that first dimension is equal to the dimensionality of the space which is one here )
    ai = hcat(a[i,:])'
    means[i,:] = kmeans(ai,2).centers
end

# normal averages for comparison
a0, Δa0, S0, _ = a_vs_central_action(h5dset, run)
Nrep = read(h5dset[run], "N_replicas")
label = L"N_s\!=\!%$Ns,N_{\!\mathrm{rep}}\!=\!%$Nrep"

plt = plot()
plot!(plt; legend = :topright, xlabel = L"u_p", ylabel = L"a_n")
LLRParsing.a_vs_central_action_plot!(plt, a0, Δa0, S0, Nt, Ns, Nrep; label,lens=false)
scatter!(plt,up,means,label="clusters")