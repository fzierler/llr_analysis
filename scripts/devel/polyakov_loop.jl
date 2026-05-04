using Pkg
Pkg.activate(".")
Pkg.instantiate()
using LLRParsing
using HDF5
using Makie
using CairoMakie
using LaTeXStrings
using PDFmerger
using ProgressMeter
using Makie.StructArrays

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end

h5file = "data_assets/sp4/all_sp4_sorted.hdf5"
h5file = "data_assets/su3/all_su3_sorted.hdf5"

h5 = h5open(h5file)
for ens in filter(!isequal("provenance"),keys(h5))    
    # lattice volume
    Nt = read(h5[ens],"Nt")
    Ns = read(h5[ens],"Ns")
    Nint = read(h5[ens],"N_replicas")
    repeats = read(h5[ens],"repeats")
    group = read(h5[ens],"group")
    # obtain fixed-a parameters
    an = read(h5["$ens/$(first(repeats))"],"an_fxa")
    S0 = read(h5["$ens/$(first(repeats))"],"S0_fxa")
    poly = read(h5["$ens/$(first(repeats))"],"poly_fxa")
    E = read(h5["$ens/$(first(repeats))"],"E_fxa")
    dS = read(h5["$ens/$(first(repeats))/Rep_0"],"dS0")
    # don't plot anything if there is no data 
    isempty(poly) && continue
    
    # for Z2 symmetric theories the polyakov loop is real
    # for ZN symmetric theories the polyakov loop is complex
    poly_re = real.(poly)
    poly_im = imag.(poly)
    poly_ang = angle.(poly)
    poly_abs = abs.(poly)
    up = E/(6Nt*Ns^3)
    up_mid = S0/(6Nt*Ns^3)

    # TODO:
    # 4) Add corresponding hightlight to shader plot
    # 6) Average over repeats, and increase dpi
    points_re = StructArray{Point2f}((vec(up), vec(poly_re)))
    points_im = StructArray{Point2f}((vec(up), vec(poly_im)))
    points3D = StructArray{Point3f}((vec(up), vec(poly_im), vec(poly_abs)))

    plotpath = "tmp_poly_plots"
    ispath(plotpath) || mkpath(plotpath)

    # select replica to highlight 
    @showprogress desc="plot Polyakov loop $ens" for rep_ind in 1:Nint

        # set up points for plotting 
        points_cplx = StructArray{Point2f}((vec(poly_re[rep_ind,:]), vec(poly_im[rep_ind,:])))

        fig = Figure(size = (600*2, 3*450))
        title = L"%$Nt\times%$(Ns)^3,~N_{\mathrm{rep}}=%$Nint"
        xlabel = L"u_p"
        set_theme!(theme_latexfonts())

        if group == "SU(3)"
            poly_label = L"\text{Im}(\ell_p)"
            ax0 = Axis(fig[1, 1]; title, xlabel = poly_label)
            ax0B = Axis(fig[1, 2]; title, xlabel = poly_label)
            ax1 = Axis(fig[2, 1]; title, xlabel, ylabel = poly_label)
            ax3 = Axis(fig[2, 2]; title, ylabel = poly_label)
            datashader!(ax1,points3D,agg = Makie.AggMean(), operation = identity, bins=2)
            datashader!(ax0B,points_cplx,colormap=[:transparent, :grey, :black], bins=2)
            hist!(ax0,vec(poly_im[rep_ind,:]))
            hist!(ax3,vec(poly_im),direction=:x)
        else
            poly_label = L"\text{Im}(\ell_p)"
            ax0 = Axis(fig[1, 1]; title, xlabel = poly_label)
            ax1 = Axis(fig[2, 1]; title, xlabel, ylabel = poly_label)
            ax3 = Axis(fig[2, 2]; title, ylabel = poly_label)
            datashader!(ax1,points_re,colormap=[:transparent, :grey, :black])
            hist!(ax0,vec(poly_re[rep_ind,:]))
            hist!(ax3,vec(poly_re),direction=:x)
        end
        ax2 = Axis(fig[3, 1]; title, xlabel, ylabel = L"a_n")
        scatter!(ax2,up_mid,an)
        vlines!(ax1,up_mid,color=:gray,alpha=0.5,linewidth=1)
        # save figure
        save(joinpath(plotpath,"$(group)_$(ens)_ind$(rep_ind).pdf"),fig)
    end
    tmp_plots = [joinpath(plotpath,"$(group)_$(ens)_ind$(rep_ind).pdf") for rep_ind in 1:Nint]
    merge_pdfs(tmp_plots, joinpath(plotpath,"$(group)_$(ens).pdf"), cleanup=true)
end