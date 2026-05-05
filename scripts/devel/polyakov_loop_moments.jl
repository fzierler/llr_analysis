using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5

h5file = "data_assets/su3/all_su3_sorted.hdf5"

h5 = h5open(h5file)
ensembles = filter(!isequal("provenance"),keys(h5))
ens = "4x20_8replicas"
an_fxa, S0, E_fxa, poly_fxa = read_fixed_data(h5,ens)

β = 5.68
repeat_id = 1
E = E_fxa[repeat_id,:,:,:]
an = an_fxa[:,repeat_id]
poly = poly_fxa[repeat_id,:,:,:]
logZ_fixed_a(E, S0, an, β)
P1 = polyakov_loop_fixed_a(E, S0, an, β, poly; f=x->abs(x)^1)
P2 = polyakov_loop_fixed_a(E, S0, an, β, poly; f=x->abs(x)^2)
P4 = polyakov_loop_fixed_a(E, S0, an, β, poly; f=x->abs(x)^4)
@show P1, P2, P4