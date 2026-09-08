using LLRParsing
using HDF5
using Peaks
using Statistics
using Plots
using LaTeXStrings
using ArgParse
using DelimitedFiles
using LsqFit
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

function apply_jackknife(obs::AbstractVector)
    N = length(obs)
    O = mean(obs)
    ΔO = sqrt(N - 1) * std(obs, corrected = false)
    return O, ΔO
end
function fix_csv_output(arr)
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
    return ßs, ß0, w0inv, w0inv_err, coeff, Cov, Nc
end

function Casmir(group, Nc)
    if (group == "SUN")
        return Nc
    elseif (group == "SPN")
        return Nc + 1
    end
end

function main(files, refw0, scale_file, out_file)
    data = readdlm(scale_file)
    mask = (data[1:end, 6] .== 0.6)
    data = data[mask, 1:end]
    Nc = (data[1:end,1])[1]
    C2 = Casmir("SPN",Nc)
    betas = (data[1:end,2])
    w0cl = (data[1:end,7])
    w0cl_err = (data[1:end,8])
    a = 1 ./ w0cl
    a_err = w0cl_err ./ (w0cl .^2)
    ß0 = 7.5;
    f(b,c) = a_beta(b,c[1],c[2],c[3],c[4],ß0,Nc,C2)
    c0 = [1.0,1.0,1.0,1.0]
    fit = curve_fit(f, betas, a,1 ./ (a_err .^ 2), c0)
    keys = ",beta,w0inv,w0inv_err,beta_0,Nc,coef,Cov"
    open(out_file, "w") do io
        println(io, keys)
        print(io,"0,")
        print(io,join(["[" * join(betas, " ") * "]"]),",")
        print(io,join(["[" * join(a, " ") * "]"]),",")
        print(io,join(["[" * join(a_err, " ") * "]"]),",")
        print(io,ß0,",")
        print(io,Nc,",")
        print(io,join(["[" * join(fit.param, " ") * "]"]),",")
        print(io,join(["[",join([join(row, " ") * " " for row in eachrow(estimate_covar(fit))]),"]"]))
    end
    return 1
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--outfile"
        help = "Where to save the output"
        required = true
        "--scalefile"
        help = "Location of the scale setting file"
        required = true
        "--refw0"
        help = "Refernce Wilson flow coefficient"
        "arg"
        help = "HDF5 files with sorted data for all Nt to be plotted"
        nargs = '+'
    end
    return parse_args(s)
end
function main()
    args = parse_commandline()
    outfile = args["outfile"]
    scalefile = args["scalefile"]
    refw0 = args["refw0"]
    files = args["arg"]
    main(files, refw0, scalefile, outfile)
    return nothing
end
main()
