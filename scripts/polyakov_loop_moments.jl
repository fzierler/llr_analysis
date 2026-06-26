using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using ArgParse
using ProgressMeter
using Statistics

function stdmean(X;dims,bin=1)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N/bin)
    return m, s
end

function main(h5file,h5file_out;nβs=100)
    h5 = h5open(h5file)
    
    runs = keys(h5)
    runs = filter(!startswith("provenance"), runs)

    for ens in runs
        if is_fixed_a_measured(h5, ens)
            an_fxa, S0, E_fxa, poly_fxa = read_fixed_data(h5,ens)

            βs = collect(range(extrema(an_fxa)...,length=nβs))
            repeats = size(an_fxa,2)

            P  = zeros(eltype(poly_fxa),length(βs),repeats)
            P1 = zeros(length(βs),repeats)
            P2 = zeros(length(βs),repeats)
            P4 = zeros(length(βs),repeats)

            @showprogress "polyakov loop: $ens" for i in eachindex(βs)
                for repeat_id in 1:repeats
                    E = E_fxa[repeat_id,:,:,:]
                    an = an_fxa[:,repeat_id]
                    poly = poly_fxa[repeat_id,:,:,:]
                    P[i,repeat_id]  = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=identity)
                    P1[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^1)
                    P2[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^2)
                    P4[i,repeat_id] = polyakov_loop_fixed_a(E, S0, an, βs[i], poly; f=x->abs(x)^4)
                end
            end
            
            lp, Δlp = stdmean(P,dims=2) 
            lp_abs, Δlp_abs = stdmean(P1,dims=2)
            lp_abs2, Δlp_abs2 = stdmean(P2,dims=2)
            lp_abs4, Δlp_abs4 = stdmean(P4,dims=2)
        else
            # if no data exists, then we save an empty array
            βs = Float64[]
            P  = ComplexF64[]
            P1 = Float64[]
            P2 = Float64[]
            P4 = Float64[]
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
        h5write(h5file_out,"$ens/lp_repeats",P)
        h5write(h5file_out,"$ens/lp_abs_repeats",P1)
        h5write(h5file_out,"$ens/lp_abs2_repeats",P2)
        h5write(h5file_out,"$ens/lp_abs4_repeats",P4)
    end
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
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    h5file_in = args["h5file_in"]
    h5file_out = args["h5file_out"]
    main(h5file_in,h5file_out)
end
main()