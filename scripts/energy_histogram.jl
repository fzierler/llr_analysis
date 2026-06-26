using LLRParsing
using HDF5
using Plots
using ArgParse
using LaTeXStrings
using PDFmerger
gr(
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 5,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
    top_margin = 0Plots.mm,
)
function read_therm_meas(h5id, name,repeat)
    N_replicas = read(h5id, "$name/N_replicas")
    llr_therm = read(h5id, "$name/$repeat/Rep_0/llr_therm")
    llr_meas = read(h5id, "$name/$repeat/Rep_0/llr_meas")
    ΔE = read(h5id, "$name/$repeat/Rep_0/dS0")
    E0 = read(h5id, "$name/$repeat/Rep_0/S0")

    E_therm = zeros(llr_therm, N_replicas)
    E_meas = zeros(llr_meas, N_replicas)
    E0 = zeros(N_replicas)

    for n in 0:(N_replicas - 1)
        E_meas[:, n + 1] = read(h5id, "$name/$repeat/Rep_$n/E_meas")[:, end]
        E_therm[:, n + 1] = read(h5id, "$name/$repeat/Rep_$n/E_therm")[:, end]
        E0[n + 1] = read(h5id, "$name/$repeat/Rep_$n/S0")[end]
    end

    # sort by central energy
    perm    = sortperm(E0)
    E_therm = E_therm[:, perm]
    E_meas  = E_meas[:, perm]

    return E_therm, E_meas, ΔE, E0
end

function main(file, name, dir; repeat = 1, extra_therm = 0 )
    h5id = h5open(file)
    repeats = read(h5id, "$name/repeats")
    Nrep = read(h5id, "$name/N_replicas")

    plt_repeat_hist = AbstractString[]
    plt_repeat_traj = AbstractString[]

    for repeat in repeats
                
        has_therm = all( [ haskey(h5id, "$name/$repeat/Rep_$n/E_therm") for n in 0:Nrep-1])
        has_meas  = all( [ haskey(h5id, "$name/$repeat/Rep_$n/E_meas")  for n in 0:Nrep-1])        

        if has_therm && has_meas
            E_therm, E_meas, ΔE, E0 = read_therm_meas(h5id, name, repeat)

            E_full = cat(E_therm,E_meas,dims=1)
            therms = size(E_therm,1) 

            nreplicas = size(E_meas, 2)
            cols = 12
            rows = cld(nreplicas,cols) 
            
            plts = [ plot(E_full[:, i], ticks = :none, label = "", title = "replica #$i") for i in axes(E_meas, 2) ]
            plts = [ vspan!(plt, [1, therms], color = :blue, alpha = 0.2, labels = "therm") for plt in plts ]
            plt1 = plot(plts..., layout = grid(cols, rows), size = (1000, 1500))
            plot!(plt1, plot_title = LLRParsing.fancy_title(name)*L"repeat $=%$repeat$")
            
            plts = [ histogram(E_meas[1+extra_therm:end, i], ticks = :none, label = "", title = "replica #$i") for i in axes(E_meas, 2) ]
            plt2 = plot(plts..., layout = grid(cols, rows), size = (1000, 1500))
            plot!(plt2, plot_title = LLRParsing.fancy_title(name)*L"repeat $=%$repeat$")

            plt_name_1 = tempname()*".pdf"
            plt_name_2 = tempname()*".pdf"
            savefig(plt1, plt_name_1)
            savefig(plt2, plt_name_2)
            push!(plt_repeat_hist,plt_name_1)
            push!(plt_repeat_traj,plt_name_2)
        end
    end
    if isempty(plt_repeat_hist) || isempty(plt_repeat_traj)
        savefig(plot(), joinpath(dir,"energy_trajectory_$name.pdf"))
        savefig(plot(), joinpath(dir,"energy_histogram_$name.pdf"))
    else
        merge_pdfs(plt_repeat_hist, joinpath(dir,"energy_trajectory_$name.pdf"), cleanup=true)
        merge_pdfs(plt_repeat_traj, joinpath(dir,"energy_histogram_$name.pdf"), cleanup=true)
    end
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--file"
        help = "HDF5 file containing the parsed data"
        required = true
        "--name"
        help = "Run to be analysed"
        required = true
        "--plot_dir"
        help = "directory in which to save the figure"
        required = true
        "--repeat"
        help = "Repeat to be plotted"
        arg_type = Int
        default = 1
        "--extra_therm"
        help = "Obtain histograms with extra thermalisation"
        arg_type = Int
        default = 0
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    file = args["file"]
    name = args["name"]
    repeat = args["repeat"]
    dir = args["plot_dir"]
    extra_therm = args["extra_therm"]
    main(file,name,dir;repeat,extra_therm)
end
main()
