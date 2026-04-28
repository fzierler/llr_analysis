using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using Plots
using LaTeXStrings
using PDFmerger
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

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end

h5file = "tmp/sp4/sp4.hdf5"
h5 = h5open(h5file)
ens = "4x20_64replicas"
rep = "0"

repeats = read(h5[ens], "repeats")
Nrep = read(h5[ens], "N_replicas")
Nt = read(h5[ens],"Nt")
Ns = read(h5[ens],"Ns")
V = Nt*Ns^3

# The following quantities are not reliably logged in the output files
# But we can deduce them from the total number of measurements
nfxa_swap = length(h5["$ens/$(first(repeats))/Rep_0/S0_fxa"])
nfxa_meas = length(h5["$ens/$(first(repeats))/Rep_0/poly"])÷nfxa_swap

S0_fxa = zeros(Nrep, nfxa_swap)
an_fxa = zeros(Nrep, nfxa_swap)
poly_fxa = zeros(ComplexF64, (Nrep, nfxa_meas, nfxa_swap))

S0_fxa_sorted = zeros(Nrep, nfxa_swap)
an_fxa_sorted = zeros(Nrep, nfxa_swap)
poly_fxa_sorted = zeros(ComplexF64, (Nrep, nfxa_meas, nfxa_swap))

for i in 1:Nrep
    S0_fxa[i, :] = read(h5["$ens/$rep/Rep_$(i - 1)"], "S0_fxa")
    an_fxa[i, :] = read(h5["$ens/$rep/Rep_$(i - 1)"], "a_fxa")
    poly_fxa[i, :, :] = read(h5["$ens/$rep/Rep_$(i - 1)"], "poly")
end

# To do: sort polyakov loop data
perm = [ sortperm(S0_fxa[:,i]) for i in axes(S0_fxa,2) ]
for (i,p) in enumerate(perm)
    S0_fxa_sorted[:,i] .= S0_fxa[p,i] 
    an_fxa_sorted[:,i] .= an_fxa[p,i] 
    poly_fxa_sorted[:,:,i] .= poly_fxa[p,:,i] 
end

up = S0_fxa_sorted[:,1]/(6V)
an = an_fxa_sorted[:,1]
poly = reshape(poly_fxa_sorted,(Nrep,nfxa_meas*nfxa_swap))

title = LLRParsing.fancy_title(ens)*", repeat #$rep"
plt_an = scatter(up,an,label="",title=title,ylabel=L"a_n",xlabel=L"u_p")

# frequency of plotting  (for smaller files)
f = 20
plt = scatter(up,real.(poly[:,1:f:end]),label="",color=:black,alpha=0.2,ms=2)
plot!(plt,ylabel=L"\ell_p",xlabel=L"u_p",title=title)

kws = (bins=:sqrt,normalize=:probability,xlabel=L"\ell_p",title=title)
plts = [histogram(real.(poly[i,1:f:end]),label="replica #$i"; kws...) for i in 1:Nrep]
for (i,p) in enumerate(plts)
    b0 = LLRParsing._highlight_replica!(deepcopy(plt),up,i; color = :green, alpha = 0.5, labels = "replica")
    a0 = LLRParsing._highlight_replica!(deepcopy(plt_an),up,i; color = :green, alpha = 0.5, labels = "replica")
    p0 = plot(p, b0, a0, layout = grid(3, 1), size = (425, 846))
end