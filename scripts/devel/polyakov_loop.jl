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
using Statistics

# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# Still need to obtain 'S' from Davids code
function logZ_fxa(S_fxa, S0, an, β)
    dS = S0[2] - S0[1]
    logρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    # figure out sorting of S0_fxa
    return VEV_exp = @. (-an[1] + β) * S_fxa[1] # +an[1]*S0[1] + logρ - log(length(S0_fxa[1])) + log(dS)
end
# For highlighting the individual energy ranges, we want to take into account 
# that the first and last interval only have a one-sided constraint
function highlight_range!(ax,xs,ind,xmin,xmax)
    if ind == 1
        x1, x2 = xmin, xs[1]
    elseif ind == length(xs)
        x1, x2 = xs[end-1], xmax
    else
        x1, x2 = xs[ind-1], xs[ind]
    end
    vspan!(ax,[x1],[x2],color=(:green, 0.5))
end
function stdmean(X;dims=1,bin=1)
    N = size(X)[dims]
    m = dropdims(mean(X;dims);dims)
    s = dropdims(std(X;dims);dims)/sqrt(N/bin)
    return m, s
end

function main(h5file)

    h5 = h5open(h5file)
    ensembles = filter(!isequal("provenance"),keys(h5))

    for ens in ensembles
        # lattice volume
        Nt = read(h5[ens],"Nt")
        Ns = read(h5[ens],"Ns")
        Nint = read(h5[ens],"N_replicas")
        Nrep = read(h5[ens],"N_repeats")
        repeats = read(h5[ens],"repeats")
        group = read(h5[ens],"group")

        # Read all central energies and coefficients a_n
        # If no fixed a-calculation was done, then an empty array will be returned 
        S0 = [ read(h5["$ens/$r"],"S0_fxa") for r in repeats]
        an = [ read(h5["$ens/$r"],"an_fxa") for r in repeats]
        # don't plot anything if there is no data 
        if any(isempty,S0) || any(isempty,an) 
            continue
        end

        # If we proceed with the plot then we have performed the measurements 
        # everywhere. Thus all entries of  S0 are idential and we can pick
        # one representative.
        # For the an's we have different values for every repeat and we rehape
        # the data into an array of size (N_intervals,N_repeats)
        S0 = first(S0)
        an = hcat(an...)

        # Determine mean and standard deviation of the mean for the 
        # fixed values of an
        an_mean, an_std = stdmean(an,dims=2)

        # Obtain number of measurements and swaps for fixed-a calculation
        s = size(h5["$ens/$(first(repeats))/E_fxa"])
        nfxa_meas, nfxa_swap = s[2:3]

        # Only these two arrays contain data that changes across repeats
        E = zeros(Nrep,Nint,nfxa_meas,nfxa_swap)
        poly = zeros(ComplexF64,(Nrep,Nint,nfxa_meas,nfxa_swap))
        for (i,r) in enumerate(repeats)
            E[i,:,:,:] .= read(h5["$ens/$r"],"E_fxa")
            poly[i,:,:,:] .= read(h5["$ens/$r"],"poly_fxa")
        end

        # for Z2 symmetric theories the polyakov loop is real
        # for ZN symmetric theories the polyakov loop is complex
        poly_re = real.(poly)
        poly_im = imag.(poly)
        poly_ang = angle.(poly)
        poly_abs = abs.(poly)
        up = E/(6Nt*Ns^3)
        up_mid = S0/(6Nt*Ns^3)
        δup = up_mid[2] - up_mid[1]

        # Number of bins to use for a histogram for a fixed energy interval 
        n_meas = size(poly,3)*size(poly,4)
        n_bins = Int(round(sqrt(n_meas)))

        # find extrema of polyakov loop for setting plot ranges
        # (enforce a symmetric intervall)
        poly_im_e = maximum(extrema(abs,poly_im))
        poly_re_e = maximum(extrema(abs,poly_re))
        poly_im_extr = (-poly_im_e,+poly_im_e)
        poly_re_extr = (-poly_re_e,+poly_re_e)

        # set up points in a enfficient data structure for the Makie datashader
        # plot. the three-dimensional data encodes is used to provide colour-
        # shading for non-Z_2 datashader plots of the Polyakov loop 
        points_re = StructArray{Point2f}((vec(up), vec(poly_re)))
        points_im = StructArray{Point2f}((vec(up), vec(poly_im)))
        points3D = StructArray{Point3f}((vec(up), vec(poly_im), vec(poly_abs)))

        plotpath = "tmp_poly_plots"
        ispath(plotpath) || mkpath(plotpath)

        # select replica to highlight 
        @showprogress desc="plot Polyakov loop $ens" for rep_ind in 1:Nint

            # set up points for plotting 
            points_cplx = StructArray{Point2f}((vec(poly_re[:,rep_ind,:,:]), vec(poly_im[:,rep_ind,:,:])))

            fig = Figure(size = (400*2, 3*350))
            title = L"%$Nt\times%$(Ns)^3,~N_{\mathrm{rep}}=%$Nint"
            xlabel = L"u_p"
            set_theme!(theme_latexfonts())

            if group == "SU(3)"
                poly_label = L"\text{Im}(\ell_p)"
                ax0 = Axis(fig[1, 1]; title, xlabel = poly_label, limits = (poly_im_extr, nothing))
                ax0B = Axis(fig[1, 2]; title, xlabel = poly_label, limits = (poly_re_extr, poly_im_extr))
                ax1 = Axis(fig[2, 1]; title, xlabel, ylabel = poly_label, limits = (extrema(up), nothing))
                ax3 = Axis(fig[2, 2]; title, ylabel = poly_label)
                datashader!(ax1,points3D,agg = Makie.AggMean(), operation = identity, binsize=3)
                datashader!(ax0B,points_cplx,colormap=[:transparent, :grey, :black], binsize=3)
                hist!(ax0,vec(poly_im[:,rep_ind,:,:]), normalization = :pdf, bins = n_bins)
                hist!(ax3,vec(poly_im),direction=:x, bins = n_bins)
            else
                poly_label = L"\text{Im}(\ell_p)"
                ax0 = Axis(fig[1, 1]; title, xlabel = poly_label, limits = (poly_re_extr, nothing))
                ax1 = Axis(fig[2, 1]; title, xlabel, ylabel = poly_label, limits = (extrema(up), nothing))
                ax3 = Axis(fig[2, 2]; title, ylabel = poly_label)
                datashader!(ax1,points_re,colormap=[:transparent, :grey, :black], binsize=3)
                hist!(ax0,vec(poly_re[:,rep_ind,:,:]), normalization = :pdf, bins = n_bins)
                hist!(ax3,vec(poly_re),direction=:x, bins = n_bins)
            end
            # Add a plot of the fixed-values of a_n averaged over repeats
            ax2 = Axis(fig[3, 1]; title, xlabel, ylabel = L"a_n", limits = (extrema(up), nothing))
            scatter!(ax2,up_mid,an_mean)
            scatter!(ax2,up_mid,an_mean)
            errorbars!(ax2, up_mid, an_mean, an_std, whiskerwidth = 10)
            # Add a histogram of the avilable values of a_n for the selected energy interval
            ax4 = Axis(fig[3, 2]; title, ylabel = L"a_n")
            hist!(ax4, an[rep_ind,:], bins = 15)
            # highlight energy interval
            highlight_range!(ax1,up_mid .+ δup/2,rep_ind,extrema(up)...)
            highlight_range!(ax2,up_mid .+ δup/2,rep_ind,extrema(up)...)
            vlines!(ax1,up_mid[1:end-1] .+ δup/2,color=:gray,alpha=0.5,linewidth=1)
            # save figure
            save(joinpath(plotpath,"$(group)_$(ens)_ind$(rep_ind).pdf"),fig)
        end
        tmp_plots = [joinpath(plotpath,"$(group)_$(ens)_ind$(rep_ind).pdf") for rep_ind in 1:Nint]
        merge_pdfs(tmp_plots, joinpath(plotpath,"$(group)_$(ens).pdf"), cleanup=true)
    end
    return nothing
end

h5file = "data_assets/sp4/all_sp4_sorted.hdf5"
h5file = "data_assets/su3/all_su3_sorted.hdf5"
main(h5file)