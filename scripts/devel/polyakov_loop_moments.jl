using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5

function main(h5file,ens;nβs=100)
    h5 = h5open(h5file)
    an_fxa, S0, E_fxa, poly_fxa = read_fixed_data(h5,ens)

    βs = range(extrema(an_fxa)...,length=nβs)
    repeats = size(an_fxa,2)

    P_resample  = zeros(eltype(poly_fxa),length(βs),repeats)
    P1_resample = zeros(length(βs),repeats)
    P2_resample = zeros(length(βs),repeats)
    P4_resample = zeros(length(βs),repeats)

    for i in eachindex(βs)
        for repeat_id in 1:repeats
            E = E_fxa[repeat_id,:,:,:]
            an = an_fxa[:,repeat_id]
            poly = poly_fxa[repeat_id,:,:,:]
            P_resample[i,repeat_id]  = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=identity)
            P1_resample[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^1)
            P2_resample[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^2)
            P4_resample[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^4)
        end
    end
    
    lp, Δlp = apply_jackknife(P_resample,dims=2) 
    lp_abs, Δlp_abs = apply_jackknife(P1_resample,dims=2)
    lp_abs2, Δlp_abs2 = apply_jackknife(P2_resample,dims=2)
    lp_abs4, Δlp_abs4 = apply_jackknife(P4_resample,dims=2)
    return lp, Δlp, lp_abs, Δlp_abs, lp_abs2, Δlp_abs2, lp_abs4, Δlp_abs4
end

h5file = "data_assets/su3/all_su3_sorted.hdf5"
ens = "4x20_8replicas"
lp, Δlp, lp_abs, Δlp_abs, lp_abs2, Δlp_abs2, lp_abs4, Δlp_abs4 = main(h5file,ens)


