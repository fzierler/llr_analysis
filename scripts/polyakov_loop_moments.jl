using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using ArgParse
using ProgressMeter

function main(h5file,h5file_out,ens;nβs=100)
    h5 = h5open(h5file)
    if ens ∈ keys(h5) && is_fixed_a_measured(h5, ens)
        an_fxa, S0, E_fxa, poly_fxa = read_fixed_data(h5,ens)

        βs = collect(range(extrema(an_fxa)...,length=nβs))
        repeats = size(an_fxa,2)

        P_resample  = zeros(eltype(poly_fxa),length(βs),repeats)
        P1_resample = zeros(length(βs),repeats)
        P2_resample = zeros(length(βs),repeats)
        P4_resample = zeros(length(βs),repeats)

        @showprogress ens for i in eachindex(βs)
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
    else
        # if no data exists, then we save an empty file
        βs = Float64[]
        P_resample  = ComplexF64[]
        P1_resample = Float64[]
        P2_resample = Float64[]
        P4_resample = Float64[]
        lp, Δlp = ComplexF64[], ComplexF64[]
        lp_abs, Δlp_abs = Float64[], Float64[]
        lp_abs2, Δlp_abs2 = Float64[], Float64[]
        lp_abs4, Δlp_abs4 = Float64[], Float64[]
    end
    h5write(h5file_out,"$ens/beta",βs)
    h5write(h5file_out,"$ens/lp",lp)
    h5write(h5file_out,"$ens/lp_abs",lp_abs)
    h5write(h5file_out,"$ens/lp_abs2",lp_abs2)
    h5write(h5file_out,"$ens/lp_abs4",lp_abs4)
    h5write(h5file_out,"$ens/Δlp",Δlp)
    h5write(h5file_out,"$ens/Δlp_abs",Δlp_abs)
    h5write(h5file_out,"$ens/Δlp_abs2",Δlp_abs2)
    h5write(h5file_out,"$ens/Δlp_abs4",Δlp_abs4)
    h5write(h5file_out,"$ens/lp_samples",P_resample)
    h5write(h5file_out,"$ens/lp_abs_samples",P1_resample)
    h5write(h5file_out,"$ens/lp_abs2_samples",P2_resample)
    h5write(h5file_out,"$ens/lp_abs4_samples",P4_resample)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file_in"
        help = "HDF5 file containing the sorted results"
        required = true
        "--h5file_out"
        help = "HDF5 file containing the polyakov loop results"
        required = true
        "--dataset"
        help = "Dataset for which the calculation will be performed"
        required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    h5file_in = args["h5file_in"]
    h5file_out = args["h5file_out"]
    ens = args["dataset"]
    main(h5file_in,h5file_out,ens)
end
main()