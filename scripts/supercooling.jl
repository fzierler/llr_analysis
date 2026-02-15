using LLRParsing
using HDF5
using ArgParse

function supercooling_all_runs(h5file, outfile, Nt)
    h5dset = h5open(h5file)
    runs = keys(h5dset)
    runs = filter(!startswith("provenance"), runs)
    if !isnothing(Nt)
        runs = filter(r -> read(h5dset[r], "Nt") == Nt, runs)
    end

    io = open(outfile, "w")
    print_provenance_csv(io)
    println(io, "group_family,Nc,n_replicas,run,T,L,t1,Δt1,t2,Δt2,tc,Δtc")
    for r in runs

        L = read(h5dset[r], "Ns")
        T = read(h5dset[r], "Nt")
        Nr = read(h5dset[r], "N_replicas")

        t1, Δt1, t2, Δt2, tc, Δtc = try
            supercooling(h5dset, r)
        catch
            NaN, NaN, NaN, NaN, NaN, NaN
        end

        println(io, "Sp,4,$Nr,$r,$T,$L,$t1,$Δt1,$t2,$Δt2,$tc,$Δtc")
    end
    close(io)
    return nothing
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--outfile"
        help = "Where to save results"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        default = 0
        arg_type = Int
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    h5file = args["h5file"]
    outfile = args["outfile"]
    Nt = args["Nt"]
    Nt = iszero(Nt) ? nothing : Nt
    return supercooling_all_runs(h5file, outfile, Nt)
end
main()
