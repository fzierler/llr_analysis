using LLRParsing
using HDF5
using DelimitedFiles

theory = "su3"
file = "data_assets/$(theory)/all_$(theory)_sorted.hdf5"
h5id = h5open(file)
runs = keys(h5id)
runs = filter(!startswith("provenance"), runs)

for r in runs
    a0, Δa0, S0, _ = a_vs_central_action(h5id, r)
    Nt = read(h5id[r], "Nt")
    Ns = read(h5id[r], "Ns")
    V = Nt * Ns * Ns * Ns
    up = S0 / (6V)

    io = open("$(theory)_$r.txt", "w")
    write(io, "beta,plaq\n")
    writedlm(io, hcat(a0, up), ',')
    close(io)
end
