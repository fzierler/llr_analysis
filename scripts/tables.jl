using LLRParsing
using HDF5
using LaTeXStrings
using ArgParse

function fmt_meas(n_meas)
    ns = unique(n_meas)
    length(ns) == 1 && return string(only(ns))
    return string(Tuple(ns))
end

function write_run_table(file, outfile)
    ispath(dirname(outfile)) || mkpath(dirname(outfile))
    fid = h5open(file)
    io = open(outfile, "w")
    print_provenance_tex(io)

    header = L"\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|} \hline
$N_t$ & $N_s$ & $u_{p}^{\rm min}$ & $u_{p}^{\rm max}$ & $N_{\rm rep}$ & $N_{\rm repeats}$ & $n_{\rm NR}$ & $n_{\rm RM}$ & $\Delta u_p$ & $n_{\rm therm}$ & $n_{\rm meas}$ \\ \hline \hline "
    footer = "\\hline \\hline
\\end{tabular}"

    runs = keys(fid)
    runs = filter(!startswith("provenance"), runs)
    println(io, header)
    for run in runs
        a0, Δa0, S0, ind = a_vs_central_action(fid, run)
        Nt = read(fid[run], "Nt")
        Ns = read(fid[run], "Ns")
        repeats = read(fid[run], "repeats")
        Nrepeats = read(fid[run], "N_repeats")
        Nreplicas = read(fid[run], "N_replicas")
        # assume that the following parameters match across repreats
        dS0 = read(fid[run], "$(first(repeats))/Rep_0/dS0")
        isrm = read(fid[run], "$(first(repeats))/Rep_0/is_rm")
        llr_meas = read(fid[run], "$(first(repeats))/Rep_0/llr_meas")
        llr_therm = read(fid[run], "$(first(repeats))/Rep_0/llr_therm")
        # count number of iterations
        up = S0 / (6Ns^3 * Nt)
        dup = dS0 / (6Ns^3 * Nt)
        nr, rm = findlast(x -> !x, isrm), length(isrm)
        nr = isnothing(nr) ? 0 : nr
        println(
            io,
            "$Nt & $Ns & $(round((minimum(up)), digits = 3)) & $(round(maximum(up), digits = 3)) & $Nreplicas & $Nrepeats & $nr & $(rm - nr) & $(round(dup, digits = 6)) & $(fmt_meas(llr_therm)) & $(fmt_meas(llr_meas)) \\\\",
        )
    end
    println(io, footer)
    return close(io)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--outfile"
        help = "Where to save the table"
        required = true
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    return write_run_table(args["h5file"], args["outfile"])
end
main()
