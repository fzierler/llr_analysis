using LLRParsing
using HDF5
using Peaks
using Statistics
using Plots
using LaTeXStrings
using ArgParse
using DelimitedFiles
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topleft,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
    palette = :Set1_5,
)
LINESTYLES = [:solid, :dot, :dash]
MARKERS = [:circle, :diamond, :dtriangle, :heptagon, :hexagon, :ltriangle, :octagon, :pentagon, :rect, :rtriangle, :star4, :star5, :star6, :star7, :star8, :utriangle ]

function a_beta(ß,c0,c1,c2,c3,ß0,Nc,C2)
    X = Nc^2 .* (1 ./ß .- 1/ß0)
    tmp1 = (ß .- ß0) .* (-12*pi^2/(11*Nc*C2))
    tmp2 = c1*X .+ c2*(X .^2) .+ c3*(X .^3)
    tmp3 = c0 .* exp.( tmp1 .+ tmp2 )
    return tmp3
end

function karsch_coeff(ß,c0,c1,c2,c3,ß0,Nc,C2)
    return 1 ./( ((-12*pi^2)/(11*Nc*C2)) .+ (-(c1*Nc^2) ./ß.^2) .+ (- (2*c2*Nc^4 .*(-(1/ß0) .+ (1 ./ß) )) ./(ß .^2)) .+ (- (3*c3*Nc^6 .*(-(1/ß0) .+ 1 ./ß) .^2) ./ (ß .^2)) )
end

function beta_latent_jackknife(fid, run; w = 5)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(fid, run)
    a_jk = LLRParsing.jackknife_resamples(a)

    beta = zeros(size(a_jk, 2))
    Pmin = zeros(size(a_jk, 2))
    Pmax = zeros(size(a_jk, 2))
    inter = zeros(size(a_jk, 2))
    Lh = zeros(size(a_jk, 2))
    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        beta[i] = LLRParsing.beta_at_equal_heights(ai, S, V)
        ups, P, ΔP = probability_density(ai, S, beta[i], V)[1:3]
        # find the two peaks of the probability distribution
        # make sure that we have found two peaks
        # find the minimum in between the peaks
        pks = findmaxima(P, w)
        @assert length(pks.indices) == 2
        i1, i2 = pks.indices
        Lh[i] = abs(ups[i1] - ups[i2])
    end
    return beta, Lh, Nt, Ns
end
function apply_jackknife(obs::AbstractVector)
    N = length(obs)
    O = mean(obs)
    ΔO = sqrt(N - 1) * std(obs, corrected = false)
    return O, ΔO
end
function fix_csv_output(arr)
    println(arr)
    return eval(Meta.parse(replace(arr, '\n' =>"")))
end

function read_scale_csv(file)
    data, header= readdlm(file,',',header = true, comments = true)
    ßs = fix_csv_output(data[1,2])
    w0inv = fix_csv_output(data[1,3])
    w0inv_err = fix_csv_output(data[1,4])
    coeff = fix_csv_output(data[1,7])
    ß0 = data[1,5]
    Nc = data[1,6]
    Cov = reshape(fix_csv_output(data[1,8]), (length(coeff), length(coeff)))
    println("Here")
    return ßs, ß0, w0inv, w0inv_err, coeff, Cov, Nc
end

function Casmir(group, Nc)
    if (group == "SUN")
        return Nc
    elseif (group == "SPN")
        return Nc + 1
    end
end

function grad_ratio(ß,ß0,c0,c1,c2,c3,Nc,C2)
    tmp1 = (((-12*pi^2)/(11*Nc*C2) - (c1*Nc^2)/ß^2 - (2*c2*Nc^4*(-(1/ß0) + 1/ß))/ß^2 - (3*c3*Nc^6*(-(1/ß0) + 1/ß)^2)/ß^2)^2*ß^2)
    return  [ 0, Nc^2/tmp1, 2*Nc^4*(1/ß-1/ß0)/tmp1,3*Nc^6*(1/ß-1/ß0)^2/tmp1]
end

function grad_ratio_ß(ß,ß0,c0,c1,c2, c3,Nc,C2)
    tmp1=( 2*c2*Nc^4)/ß^4
    + (6*c3*Nc^6*(-(1/ß0) + 1/ß))/ß^4
    + (2*c1*Nc^2)/ß^3
    + (4*c2*Nc^4*(-(1/ß0) + 1/ß))/ß^3
    + (6*c3*Nc^6*(-(1/ß0) + 1/ß)^2)/ß^3
    tmp2 = ((-12*pi^2)/(11*Nc*C2)
            -(c1*Nc^2)/ß^2
            -(2*c2*Nc^4*(-(1/ß0) + 1/ß))/ß^2
            -(3*c3*Nc^6*(-(1/ß0) + 1/ß)^2)/ß^2 )^2
    return -tmp1/tmp2
end

function latent_heat(up,Δup,Nt,ß,Δß,Nc,C2,ß0,c0,c1,c2,c3,Cov)
    k = karsch_coeff(ß,c0,c1,c2,c3,ß0,Nc,C2)
    Lh = -6 * Nt^4 * up * k
    println(Lh)
    gradr = grad_ratio(ß,ß0,c0,c1,c2,c3,Nc,C2)
    gradß = grad_ratio_ß(ß,ß0,c0,c1,c2, c3,Nc,C2)
    ΔLh = Lh * sqrt((Δup/up)^2 + ( transpose(gradr) * Cov * gradr )/ k^2 + (Δß * gradß )^2  / k^2)
    return Lh, ΔLh
end

function main(files, scale_file, plt_name)
    plt = plot(; ylabel = L"\langle \Delta u_p \rangle_{\beta_{CV}}", xlabel = L"(N_t/N_s)^3")
    Nt = 0
    ßs, ß0, w0inv, w0inv_err, coeff, Cov, Nc = read_scale_csv(scale_file)
    c0, c1, c2, c3 =  coeff[1], coeff[2], coeff[3], coeff[4]
    ß0 = 7.5
    C2 = Casmir("SPN", Nc)
    for (j,file) in enumerate(files)
        fid = h5open(file)
        runs = keys(fid)
        runs = filter(!startswith("provenance"), runs)
        x, ß, Δß, up, Δup, Lh, ΔLh  = zeros(length(runs)), zeros(length(runs)), zeros(length(runs)), zeros(length(runs)), zeros(length(runs)), zeros(length(runs)), zeros(length(runs))
        for (i, r) in enumerate(runs)
            beta, latent, Nt, Ns = beta_latent_jackknife(fid, r)
            ß[i], Δß[i] = apply_jackknife(beta)
            up[i], Δup[i] = apply_jackknife(latent)
            Lh[i],ΔLh[i]  = latent_heat(up[i], Δup[i],Nt,ß[i], Δß[i],Nc,C2,ß0,c0,c1,c2,c3,Cov)
            x[i] = inv(Ns / Nt)
        end
        plot!(plt, x .^ 3, Lh, yerr = ΔLh, markershape = MARKERS[j], markeralpha = 0.7, label = L"N_t=%$Nt")
    end
    plot!(plt; ylims = (minimum(ylims(plt)), maximum(ylims(plt))))
    plot!(plt; xlims = (0, maximum(xlims(plt))))
    savefig(plt,  plt_name)
    return plt
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--plotfile"
        help = "Where to save the plot"
        required = true
        "--scalefile"
        help = "Location of the scale setting file"
        required = true
        "arg"
        help = "HDF5 files with sorted data for all Nt to be plotted"
        required = true
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    plotfile = args["plotfile"]
    scalefile = args["scalefile"]
    files = args["arg"]
    main(files, scalefile, plotfile)
    return nothing
end
main()
