using LLRParsing
using HDF5

h5file = "data_assets/SU3_Nt4_sorted.hdf5"
h5file = "tmp/SU3_Nt4.hdf5"
h5   = h5open(h5file)

S0 = [ last(read(h5["4x20_8replicas/0/Rep_$i"],"S0")) for i in 0:7]
an = [ last(read(h5["4x20_8replicas/0/Rep_$i"],"a")) for i in 0:7]
S0_fxa = [ read(h5["4x20_8replicas/0/Rep_$i"],"S0_fxa") for i in 0:7]
an_fxa = [ read(h5["4x20_8replicas/0/Rep_$i"],"a_fxa") for i in 0:7]
poly_fxa = [ read(h5["4x20_8replicas/0/Rep_$i"],"poly") for i in 0:7]

S0_fxa[1]
an_fxa[1]
poly_fxa[1]

perm = sortperm(S0)
permute!(an,perm)
permute!(S0,perm)
permute!(poly_fxa,perm)

function logZ_fxa(S0_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1],S0,dS,an)

    # figure out sorting of S0_fxa
    VEV_exp = @. (-an[1] + β)*S0_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
    @show S0_fxa[:][1]
    @show VEV_exp[1]
end

β = 5.68
logZ_fxa(S0_fxa, S0, an, β)

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code

# The following quantities are not reliably logged in the output files
# llr:sfreq_fxa = 20
# llr:nfxa = 2000
#
# Note that size of the arrays parsed from a single replica
# length(an[1]) = 2000
# length(S0[1]) = 2000
# length(poly[1]) = 40000
#
