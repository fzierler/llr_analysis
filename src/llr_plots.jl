function a_vs_central_action_repeats(h5dset, run; ind = nothing)
    N_replicas = read(h5dset[run], "N_replicas")
    N_repeats = read(h5dset[run], "N_repeats")
    a::Matrix{Float64} = zeros(N_replicas, N_repeats)
    S::Matrix{Float64} = zeros(N_replicas, N_repeats)
    repeat_indices = Int[]
    # read all last elements for a and the central action
    for i in 1:N_replicas, j in 1:N_repeats
        if haskey(h5dset[run], "$(j - 1)")
            n_steps = length(h5dset[run]["$(j - 1)/Rep_$(i - 1)/is_rm"])
            if n_steps > 0
                ind = isnothing(ind) ? n_steps : min(ind, n_steps)
                ind = max(1, ind)
                a[i, j] = h5dset[run]["$(j - 1)/Rep_$(i - 1)/a_sorted"][ind]
                S[i, j] = h5dset[run]["$(j - 1)/Rep_$(i - 1)/S0_sorted"][ind]
                append!(repeat_indices, j)
            end
        end
    end
    repeat_indices::Vector{Int} = unique(repeat_indices)
    return a[:, repeat_indices], S[:, repeat_indices], ind
end
function a_vs_central_action(h5dset, run; ind = nothing)
    a, S, ind = a_vs_central_action_repeats(h5dset, run; ind)
    N = size(a)[2]
    S0 = S[:, 1]
    a0 = dropdims(mean(a, dims = 2), dims = 2)
    Δa0 = dropdims(std(a, dims = 2), dims = 2) / sqrt(N)
    return a0, Δa0, S0, ind
end
function a_vs_central_action_plot(h5dset, runs; indices)
    plt = plot()
    for run in runs
        for ind in indices
            a0, Δa0, S0, ind = a_vs_central_action(h5dset, run; ind)
            Nt = read(fid[run], "Nt")
            Ns = read(fid[run], "Ns")
            V = Ns^3 * Nt
            up = S0 / (6V)
            plot!(
                plt,
                up,
                a0,
                yerr = Δa0,
                marker = :auto,
                label = "$(Nt)x$(Ns): ΔE=$(round(2(S0[2] - S0[1]) / 6V, sigdigits = 1)) (Nr+RM steps=$ind)",
            )
        end
    end
    return plt
end
function _lens_location(a0, Δa0, up)
    # identify approximate midpoint by the maximum of the variance of a_n
    # then estimate peak and dip position by looking for
    #   1) the maximum before the midpoint
    #   2) the minimum after the midpoint
    mid, m_ind = findmax(Δa0)
    peak, p_ind = findmax(a0[1:(m_ind + 1)])
    dip, d_ind = findmin(a0[m_ind:end])
    delta_peak_dip = max((m_ind - p_ind), d_ind)
    ylms = [dip, peak] .+ [dip - peak, peak - dip] ./ 3
    xinds = [p_ind - 3delta_peak_dip ÷ 2, m_ind + 7delta_peak_dip ÷ 2]
    xlms = up[xinds]
    return ylms, xlms, [dip, peak], up[[p_ind, m_ind + d_ind]]
end
function a_vs_central_action_plot!(plt, h5dset, run; index = nothing, kws...)
    a0, Δa0, S0, ind = a_vs_central_action(h5dset, run; ind = index)
    Nt = read(h5dset[run], "Nt")
    Ns = read(h5dset[run], "Ns")
    Nrep = read(h5dset[run], "N_replicas")
    return a_vs_central_action_plot!(plt, a0, Δa0, S0, Nt, Ns, Nrep; kws...)
end
function a_vs_central_action_plot!(
        plt,
        a0,
        Δa0,
        S0,
        Nt,
        Ns,
        replicas;
        label = LLRParsing.fancy_title(Nt, Ns, replicas),
        highlight_index = nothing,
        lens = true,
    )
    V = Ns^3 * Nt
    up = S0 / (6V)
    plot!(plt, up, a0, yerr = Δa0, marker = :auto; label)
    if !isnothing(highlight_index)
        mid = up[highlight_index]
        del = (up[highlight_index + 1] - up[highlight_index]) / 2
        vspan!(plt, [mid - del, mid + del], color = :green, alpha = 0.7, labels = "replica")
    end
    if lens
        ylms, xlms, yticks, xticks = _lens_location(a0, Δa0, up)
        xticks = round.(xticks, sigdigits = 4)
        yticks = round.(yticks, sigdigits = 5)
        lens!(plt, xlms, ylms, inset = (1, bbox(0.6, 0.5, 0.38, 0.38)))
        plot!(; subplot = 2, yticks, xticks)
    end
    return plt
end
function a_variance_vs_central_action_plot!(
        plt,
        h5dset,
        run;
        index = nothing,
        highlight_index = nothing,
    )
    a0, Δa0, S0, ind = a_vs_central_action(h5dset, run; ind = index)
    Nt = read(h5dset[run], "Nt")
    Ns = read(h5dset[run], "Ns")
    V = Ns^3 * Nt
    up = S0 / (6V)
    plot!(plt, up, zero(a0), ribbon = (zero(Δa0), Δa0), label = L"\Delta a_n")
    return if !isnothing(highlight_index)
        mid = up[highlight_index]
        del = (up[highlight_index + 1] - up[highlight_index]) / 2
        vspan!(plt, [mid - del, mid + del], color = :green, alpha = 0.8, labels = "replica")
    end
end
function a_trajectory(h5dset, run; replica = 0)
    N_replicas = read(h5dset[run], "N_replicas")
    repeats = read(h5dset[run], "repeats")
    @assert replica < N_replicas
    a = [Float64[] for _ in repeats]
    repeat_indices = Int[]
    for (j, repeat) in enumerate(repeats)
        if haskey(h5dset[run], repeat)
            a[j] = vec(h5dset[run]["$repeat/Rep_$(replica)/a_sorted"][])
            append!(repeat_indices, j)
        end
    end
    return a[repeat_indices]
end
function plot_a_trajectory_repeat!(plt, h5dset, run, repeat, replica)
    a = a_trajectory(h5dset, run; replica)
    plot!(plt, a[repeat], label = "", lw = 1, markersize = 1, marker = :circle, msw = 0.1)
    return plt
end
function plot_a_trajectory_all!(plt, h5dset, run, replica)
    a = a_trajectory(h5dset, run; replica)
    plot!(plt, a, label = "", lw = 1, markersize = 1, marker = :circle, msw = 0.1)
    return plt
end
function plot_nr_rm_shading!(plt, h5dset, run, repeat, replica)
    isrm = read(h5dset[run], "$repeat/Rep_$replica/is_rm")
    nr, rm = findlast(x -> !x, isrm), length(isrm)
    nr = isnothing(nr) ? 0 : nr - 1
    vspan!(plt, [1, nr + 1], color = :green, alpha = 0.2, labels = "NR")
    vspan!(plt, [nr + 1, rm], color = :blue, alpha = 0.2, labels = "RM")
    return plt
end
function plot_a_repeat_average!(plt, h5dset, run; replica)
    a = a_trajectory(h5dset, run; replica)
    # remove repeats that have a zero-length trajectory
    filter!(x -> length(x) > 0, a)
    a_last = last.(a)
    N_repeats = length(a_last)
    a0, Δa = mean(a_last), std(a_last) / sqrt(N_repeats)
    hspan!(plt, [a0 - Δa, a0 + Δa], color = :black, alpha = 0.8, labels = L"a_n")
    return plt
end
fancy_title(Lt, Ls, Nint) = L"%$Lt\!\times\!%$(Ls)^3\!, N_{\!\mathrm{rep}}\!\!=%$Nint"
function fancy_title(run)
    rx = r"([0-9]+)x([0-9]+)_([0-9]+)replicas"
    m = match(rx, run).captures
    Lt, Ls, Nint = m
    return fancy_title(Lt, Ls, Nint)
end
function full_trajectory_plot(h5dset, run, repeat_id, replica_id; lens = true)
    plt1 = plot_a_trajectory_repeat!(plot(), h5dset, run, repeat_id, replica_id)
    plt2 = plot_a_trajectory_all!(plot(), h5dset, run, replica_id)
    plt3 = a_vs_central_action_plot!(
        plot(),
        h5dset,
        run;
        lens,
        index = nothing,
        highlight_index = replica_id,
    )
    plt4 = a_variance_vs_central_action_plot!(
        plot(),
        h5dset,
        run;
        index = nothing,
        highlight_index = replica_id,
    )
    for plt in [plt1, plt2]
        plot_nr_rm_shading!(plt, h5dset, run, repeat_id, replica_id)
        plot_a_repeat_average!(plt, h5dset, run; replica = replica_id)
    end
    plot!(plt1, ylims = ylims(plt2), title = fancy_title(run), ylabel = L"a_n^{(m)}")
    plot!(
        plt2,
        ylims = ylims(plt2),
        title = "",
        xlabel = "NR/RM iteration m",
        ylabel = L"a_n^{(m)}",
    )
    plot!(plt3, title = "", xlabel = "", ylabel = L"a_n")
    plot!(
        plt4,
        title = "",
        xlabel = L"\textrm{central}~\textrm{plaquette}~u_p",
        ylabel = L"\Delta a_n",
    )
    if lens
        plot!(plt3, subplot = 2, left_margin = 0Plots.mm, tickfontsize = 8, ylabel = "")
    end
    l = grid(4, 1, heights = [0.29, 0.29, 0.29, 0.13])
    plt = plot(plt1, plt2, plt3, plt4, layout = l, size = (500, 1000))
    return plt
end
