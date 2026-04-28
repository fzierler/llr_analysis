using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end

h5file = "tmp/su3/su3.hdf5"
h5 = h5open(h5file)
ens = "4x20_8replicas"
Nrep = read(h5[ens], "N_replicas")

# The following quantities are not reliably logged in the output files
# llr:sfreq_fxa = 20
# llr:nfxa = 2000
nfxa_meas = 20
nfxa_swap = 2000

S0_fxa = zeros(Nrep, nfxa_swap)
an_fxa = zeros(Nrep, nfxa_swap)
poly_fxa = zeros(ComplexF64, (Nrep, nfxa_meas, nfxa_swap))

for i in 1:Nrep
    S0_fxa[i, :] = read(h5["$ens/0/Rep_$(i - 1)"], "S0_fxa")
    an_fxa[i, :] = read(h5["$ens/0/Rep_$(i - 1)"], "a_fxa")
    poly_fxa[i, :, :] = read(h5["$ens/0/Rep_$(i - 1)"], "poly")
end
