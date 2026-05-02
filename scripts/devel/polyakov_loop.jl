using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using Plots
using LaTeXStrings
using PDFmerger
using ProgressMeter
gr(
    size = (425, 282),
    fontfamily = "Computer Modern",
    legend = :topright,
    frame = :box,
    titlefontsize = 10,
    legendfontsize = 7,
    tickfontsize = 7,
    labelfontsize = 10,
    left_margin = 0Plots.mm,
)

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end
function plot_poly_overview(an,up,poly,plotpath,plotname,repeat,ens;f1,f2)
    Nrep = length(an)

    title = LLRParsing.fancy_title(ens)*", repeat #$repeat"
    plt_an = scatter(up,an,label="",title=title,ylabel=L"a_n",xlabel=L"u_p")

    plt = scatter(up,poly[:,1:f2:end],label="",color=:black,alpha=0.2,ms=2)
    plot!(plt,ylabel=L"\ell_p",xlabel=L"u_p",title=title)

    kws = (bins=:sqrt,normalize=:probability,xlabel=L"\ell_p",title=title)
    plts = [histogram(poly[i,1:f1:end],label="replica #$i"; kws...) for i in 1:Nrep]

    tmpfiles = [ tempname()*"_$i.pdf" for i in eachindex(plts) ]
    @showprogress desc=ens for (i,p) in enumerate(plts)
        b0 = LLRParsing._highlight_replica!(deepcopy(plt),up,i; color = :green, alpha = 0.5, labels = "replica")
        a0 = LLRParsing._highlight_replica!(deepcopy(plt_an),up,i; color = :green, alpha = 0.5, labels = "replica")
        p0 = plot(p, b0, a0, layout = grid(3, 1), size = (425, 846))
        savefig(p0,tmpfiles[i])
    end
    # combine these plots into a single file
    merge_pdfs(tmpfiles,joinpath(plotpath,plotname),cleanup=true)
end

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
    up = S0/(6Nt*Ns^3)

    # for Z2 symmetric theories the polyakov loop is real
    # for ZN symmetric theories the polyakov loop is complex
    poly = real.(poly)

    # frequency of plotting  (for smaller files)
    f1 = 50
    f2 = 500
    plotpath = "plots_hist_new"
    plotname = "polyakov_loop_su3_$(ens)_scatter.pdf"
    ispath(plotpath) || mkpath(plotpath)
    plot_poly_overview(an,up,poly,plotpath,plotname,first(repeats),ens;f1,f2)
end