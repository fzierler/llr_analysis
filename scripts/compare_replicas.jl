using LLRParsing
using LaTeXStrings
using HDF5
using Plots
using ArgParse
using Peaks
gr(
    size = (350, 300),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 9,
    tickfontsize = 9,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

a_vs_central_action_plot(h5id, runs::Vector, Nt, Ns; kws...) =
    a_vs_central_action_plot!(plot(), h5id, runs, Nt, Ns; kws...)
function a_vs_central_action_plot!(plt, h5id, runs::Vector, Nt, Ns; kws...)
    # default plot limits
    xmin, xmax = +Inf, -Inf
    ymin, ymax = +Inf, -Inf
    for run in runs
        Nt0 = read(h5id[run], "Nt")
        Ns0 = read(h5id[run], "Ns")
        if Nt0 == Nt && Ns0 == Ns
            a0, Δa0, S0, _ = a_vs_central_action(h5id, run)
            Nrep = read(h5id[run], "N_replicas")
            up = @. S0 / (6 * Ns^3 * Nt)
            label = L"N_{\!\mathrm{rep}}\!=\!%$Nrep"
            LLRParsing.a_vs_central_action_plot!(
                plt,
                a0,
                Δa0,
                S0,
                Nt,
                Ns,
                Nrep;
                label,
                kws...,
            )
            # find useful plot limits for the volume comparison
            p_ind = findmaxima(a0, 5).indices
            m_ind = findminima(a0, 5).indices
            if length(m_ind) == length(p_ind) == 1
                p_ind = only(p_ind)
                m_ind = only(m_ind)
                δ = m_ind - p_ind
                xmin, xmax = min(xmin, up[p_ind - δ]), max(xmax, up[m_ind + 2δ])
                ymin, ymax = min(ymin, a0[p_ind - δ]), max(ymax, a0[m_ind + 2δ])
            else
                ind = findmax(Δa0)[2]
                δind = min(ind, length(a0) - ind) ÷ 4
                xmin, xmax = min(xmin, up[ind - δind]), max(xmax, up[ind + δind])
                ymin, ymax = min(ymin, a0[ind - δind]), max(ymax, a0[ind + δind])
            end
        end
    end
    plot!(
        plt,
        title = L"N_t\times N_s^3\!=\!%$Nt\times %$Ns^3",
        xlims = (xmin, xmax),
        ylims = (ymin, ymax),
    )
    return plt
end
function an_action_volumes_replicas(file, plotdest, Nt, Ns)
    ispath(dirname(plotdest)) || mkpath(dirname(plotdest))
    h5id = h5open(file)
    runs = keys(h5id)
    runs = filter(!startswith("provenance"), runs)
    plt = a_vs_central_action_plot(h5id, runs, Nt, Ns, lens = false)
    plot!(plt; legend = :bottomright, xlabel = L"u_p", ylabel = L"a_n")
    return savefig(plt, plotdest)
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--h5file"
        help = "HDF5 file containing the sorted results"
        required = true
        "--plot_file"
        help = "Where to save the plot"
        required = true
        "--Nt"
        help = "Nt of the runs to be plotted of the plot"
        required = true
        arg_type = Int
        "--Ns"
        help = "Ns of the runs to be plotted of the plot"
        required = true
        arg_type = Int
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["h5file"]
    plotdst = args["plot_file"]
    Nt = args["Nt"]
    Ns = args["Ns"]
    return an_action_volumes_replicas(file, plotdst, Nt, Ns)
end
main()
