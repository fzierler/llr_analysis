using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using Makie
using CairoMakie
using LaTeXStrings
using PDFmerger
using ProgressMeter
using Makie.StructArrays

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end

h5file = "data_assets/su3/all_su3_sorted.hdf5"
h5file = "data_assets/sp4/all_sp4_sorted.hdf5"
h5 = h5open(h5file)
for ens in filter(!isequal("provenance"),keys(h5))    
    # lattice volume
    Nt = read(h5[ens],"Nt")
    Ns = read(h5[ens],"Ns")
    repeats = read(h5[ens],"repeats")
    # obtain fixed-a parameters
    an = read(h5["$ens/$(first(repeats))"],"an_fxa")
    S0 = read(h5["$ens/$(first(repeats))"],"S0_fxa")
    poly = read(h5["$ens/$(first(repeats))"],"poly_fxa")
    E = read(h5["$ens/$(first(repeats))"],"E_fxa")
    # don't plot anything if there is no data 
    isempty(poly) && continue
    
    # for Z2 symmetric theories the polyakov loop is real
    # for ZN symmetric theories the polyakov loop is complex
    poly_re = real.(poly)
    poly_im = imag.(poly)
    poly_ang = angle.(poly)
    poly_abs = abs.(poly)
    up = E/(6Nt*Ns^3)

    # set up points for plotting 
    points = StructArray{Point2f}((vec(up), vec(poly_re)))

    fig = Figure()
    set_theme!(theme_latexfonts())
    ax = Axis(fig[1, 1], title = "Title")
    datashader!(ax,points,colormap=[:transparent, :grey, :black])
    save("$ens.pdf",fig)
end