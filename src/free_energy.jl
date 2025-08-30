function thermodynamic_potentials_repeats(h5dset, run; kws...)
    a, S0_all, ind = LLRParsing.a_vs_central_action_repeats(h5dset, run)
    Ns = read(h5dset[run], "Ns")
    Nt = read(h5dset[run], "Nt")
    V = Nt * Ns^3
    S0 = dropdims(unique(S0_all, dims = 2), dims = 2)
    return thermodynamic_potentials_repeats(a, S0, V; kws...)
end
function thermodynamic_potentials(h5dset, run; kws...)
    a, S0_all, ind = LLRParsing.a_vs_central_action_repeats(h5dset, run)
    Ns = read(h5dset[run], "Ns")
    Nt = read(h5dset[run], "Nt")
    V = Nt * Ns^3
    S0 = dropdims(unique(S0_all, dims = 2), dims = 2)
    return thermodynamic_potentials(a, S0, V; kws...)
end
function thermodynamic_potentials_repeats(a, S0, V; s0 = 0.0)
    @. S0 = S0 / V
    dS = S0[2] - S0[1]
    E = copy(S0)
    s, t, f, u = zero(a), zero(a), zero(a), zero(a)
    replicas = size(a, 1)
    repeats = size(a, 2)
    for r in 1:repeats
        a0 = a[:, r]
        cs = cumsum(a0)
        for i in 1:replicas
            logρ = log_rho(E[i], S0, dS, a0; cumsum_a = cs)
            s[i, r] = logρ
            t[i, r] = 1 / a0[i]
        end
        # Perform additive constant to the entropy
        min_s = minimum(s[:, r])
        @. s[:, r] = s[:, r] - min_s
        for i in 1:replicas
            u = 6 - E[i]
            f[i, r] = u - t[i, r] * (s[i, r] + s0)
        end
    end
    return t, f, s
end
function thermodynamic_potentials(a, S0, V; kws...)
    t_r, f_r, s_r = thermodynamic_potentials_repeats(a, S0, V; kws...)
    N = size(t_r, 2)
    t = dropdims(mean(t_r, dims = 2), dims = 2)
    f = dropdims(mean(f_r, dims = 2), dims = 2)
    s = dropdims(mean(s_r, dims = 2), dims = 2)
    Δt = dropdims(std(t_r, dims = 2), dims = 2) ./ sqrt(N)
    Δf = dropdims(std(f_r, dims = 2), dims = 2) ./ sqrt(N)
    Δs = dropdims(std(s_r, dims = 2), dims = 2) ./ sqrt(N)
    return t, Δt, f, Δf, s, Δs
end
function plot_free_energies(file, plotdir, critical_entropy)
    h5dset = h5open(file)
    runs = keys(h5dset)
    close(h5dset)
    for r in runs
        plot_free_energy(file, joinpath(plotdir, "$r.pdf"), r, critical_entropy)
    end
    return
end
function plot_free_energy(file, plotfile, run, critical_entropy)
    h5dset = h5open(file)

    a, Δa, S0, _ = a_vs_central_action(h5dset, run)
    t, Δt, f, Δf, s, Δs = thermodynamic_potentials(h5dset, run)
    # show slope of the interface region
    try 
        pks = only(findmaxima(a, 5).indices)
        mns = only(findminima(a, 5).indices)
        t1, t2 = t[pks], t[mns]
        βc = LLRParsing.beta_at_equal_heights(h5dset, run)
        ind_tc = argmin(i -> abs(a[i] - βc), pks:mns)
        
        # add constant s0 such that the entropy matches at criticality
        t, Δt, f, Δf, s, Δs = thermodynamic_potentials(h5dset, run; s0 = critical_entropy - s[ind_tc])
        
        # find crossing in free energy
        r1 = 1:pks
        r2 = mns:length(a)
        perm1 = sortperm(t[r1])
        perm2 = sortperm(t[r2])
        itp1 = Interpolator(t[r1][perm1], f[r1][perm1])
        itp2 = Interpolator(t[r2][perm2], f[r2][perm2])
        
        g(x) = itp1(x) - itp2(x)
        t1, t2 = extrema(filter(t -> isfinite(g(t)), vcat(t[r1], t[r2])))
        tc = find_zero(g, (t1, t2))
        fc = itp1(tc)
        
        # set plot limits
        tmin, tmax = extrema((t[pks], t[mns]))
        δt = tmax - tmin
        tmin, tmax = tmin - δt / 4, tmax + δt / 4
        fmin, fmax = extrema((f[pks], f[mns]))
        δf = fmax - fmin
        fmin, fmax = extrema((itp2(tmin), itp1(tmax), fmin - δf / 3, fmax + δf / 3))
        
        # rescale free energy for nicer, centred plots
        scale = 10^6
        @. f = f - fc
        @. f = f * scale
        @. Δf = Δf * scale
        fmin, fmax = (fmin - fc) * scale, (fmax - fc) * scale
        
        ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
        plt = plot(title = LLRParsing.fancy_title(run))
        plot!(; ylabel = L"(f - f_c^+ )/ 10^{-6}", xlabel = L"t = 1/a_n")
        plot!(plt, t, f, xerr = Δt, yerr = Δf, ms = 1, label = "")
        plot!(plt, ylims = (fmin, fmax), xlims = (tmin, tmax), xformatter = :plain)
        savefig(plt, plotfile)
    catch 
        ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
        plt = plot(title = LLRParsing.fancy_title(run))
        plot!(; ylabel = L"(f - f_c^+ )/ 10^{-6}", xlabel = L"t = 1/a_n")
        savefig(plt, plotfile)
    end
    return close(h5dset)
end
function plot_entropy(file, plotfile, critical_entropy)
    h5dset = h5open(file)
    runs = keys(h5dset)
    runs = filter(!startswith("provenance"), runs)
    Nt = read(h5dset, joinpath(first(runs), "Nt"))
    plt = plot(title = L"$N_t = %$Nt$")
    for r in runs
        a, Δa, S0, _ = a_vs_central_action(h5dset, r)
        t, Δt, f, Δf, s, Δs = thermodynamic_potentials(h5dset, r)

        # show slope of the interface region
        pks = findmaxima(a, 5).indices
        mns = findminima(a, 5).indices

        if length(pks)==1 && length(mns)==1
            pks = only(pks)
            mns = only(mns)

            try
                βc = LLRParsing.beta_at_equal_heights(h5dset, r)
                ind_tc = argmin(i -> abs(a[i] - βc), pks:mns)
                s0 = critical_entropy - s[ind_tc]
                @. s = s + s0

                plot!(plt, t, s, xerr = Δt, yerr = Δs, ms = 1, label = LLRParsing.fancy_title(r))
            catch
            end
        end
    end
    ispath(dirname(plotfile)) || mkpath(dirname(plotfile))
    plot!(plt, legend = :bottomright, xlabel = L"t = 1/a_n", ylabel = L"$s$")
    return savefig(plt, plotfile)
end
