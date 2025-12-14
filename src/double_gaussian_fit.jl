Gaussian(x, A, μ, σ) = A * exp(-(x - μ)^2 / σ^2 / 2)
DoubleGaussian(x, A1, μ1, σ1, A2, μ2, σ2) =
    Gaussian(x, A1, μ1, σ1) + Gaussian(x, A2, μ2, σ2)
@. modelDG(x, p) = DoubleGaussian(x, p[1], p[2], p[3], p[4], p[5], p[6])
function initial_param_double_gaussian(ups, P; w = 5)
    δups = ups[2] - ups[1]
    pks = findmaxima(P, w)
    pks = peakproms(pks)
    pks = peakwidths(pks)
    p0 = [
        pks.heights[1],
        ups[pks.indices[1]],
        pks.widths[1] * δups,
        pks.heights[2],
        ups[pks.indices[2]],
        pks.widths[2] * δups,
    ]
    return p0
end
function fitting_range_double_gaussian(P; w = 5)
    # determine fitting range:
    # Include half the points between peak and minimum
    pks = findmaxima(P, w)
    mns = findminima(P, w)
    # After discussion: Only fit up to the peaks
    δ = 0
    range1 = 1:(pks.indices[1] + δ)
    range2 = (pks.indices[2] - δ):length(P)
    range = vcat(collect(range1), collect(range2))
    return range
end
function fit_double_gaussian(
        ups::Vector{T},
        P::Vector{T},
        covP::Matrix{T};
        kws...,
    ) where {T}
    p0 = initial_param_double_gaussian(ups, P; kws...)
    r = fitting_range_double_gaussian(P; kws...)
    # Remove negative eigenvalues of the correlation Matrix by hand
    # Thse are due to numerics, since we know that the correlation
    # matrix is at least semidefinite.
    λ = eigmin(Hermitian(covP[r, r]))
    fit = curve_fit(modelDG, ups[r], P[r], p0)
    return fit
end
function fit_double_gaussian(fid::HDF5.File, run, beta; kws...)
    ups, P, ΔP, covP, V, dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups, P, covP; kws...)
    return fit
end
function jackknife_resamples(obs)
    nops, nconf = size(obs)
    samples = similar(obs)
    # Strategy: Sum over all configs then remove one
    tmp = dropdims(sum(obs, dims = 2), dims = 2)
    for index in 1:nconf
        @. samples[:, index] = (tmp - obs[:, index]) / (nconf - 1)
    end
    return samples
end
function histogram_jackknife_fit(fid, run)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(fid, run)
    # estimate the covariance matrix of the probability distribution by calculating
    # it from the ensemble averages. I will use it as the weight in the least squares fit.
    a_jk = jackknife_resamples(a)
    fitted = similar(a_jk)
    ups = similar(S)
    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        beta = LLRParsing.beta_at_equal_heights(ai, S, V)
        # Note, the standard deviation and covariance matrix returned are useless
        # because we are only considering one resample at a time
        ups, P, _, covP, V, dS = LLRParsing.probability_density(ai, S, beta, V)
        fit = LLRParsing.fit_double_gaussian(ups, P, covP)
        data = 6 .* V .* LLRParsing.modelDG(ups, fit.param)
        fitted[:, i] = data
    end
    f, Δf = apply_jackknife(fitted, dims = 2)
    return ups, f, Δf
end
function βc_jackknife(fid, run; kws...)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(fid, run)
    a_jk = jackknife_resamples(a)
    beta = zeros(size(a_jk)[2])
    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        beta[i] = LLRParsing.beta_at_equal_heights(ai, S, V; kws...)
    end
    βc, Δβc = apply_jackknife(beta)
    return βc, Δβc
end
function histogram_jackknife_fit(fid, run, beta)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(fid, run)
    # estimate the covariance matrix of the probability distribution by calculating
    # it from the ensemble averages. I will use it as the weight in the least squares fit.
    covP = LLRParsing.probability_density(a, S, beta, V)[4]
    a_jk = jackknife_resamples(a)
    fitted = similar(a_jk)
    ups = similar(S)
    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        # Note, the standard deviation and covariance matrix returned are useless
        # because we are only considering one resample at a time
        ups, P, _, _, V, dS = LLRParsing.probability_density(ai, S, beta, V)
        fit = LLRParsing.fit_double_gaussian(ups, P, covP)
        data = 6 .* V .* LLRParsing.modelDG(ups, fit.param)
        fitted[:, i] = data
    end
    f, Δf = apply_jackknife(fitted, dims = 2)
    return ups, f, Δf
end
function plot_double_gaussian_fit!(plt, fid, run, beta; kws...)
    ups, P, ΔP, covP, V, dS = probability_density(fid, run, beta)
    fit = fit_double_gaussian(ups, P, covP; kws...)
    fitted = 6 .* V .* modelDG(ups, fit.param)
    plot!(plt, ups, fitted, label = "double Gaussian fit", lw = 3)
    return plt
end
function plot_double_gaussian_fit_difference(plt, fid, run; kws...)
    ups, f, Δf = histogram_jackknife_fit(fid, run)
    βc = beta_at_equal_heights(fid, run; A1 = 1, A2 = 1)
    ups, P, ΔP, covP, V, dS = probability_density(fid, run, βc)
    d, Δd = P * 6V, ΔP * 6V
    return plot!(plt, ups, d - f; ribbon = sqrt.(Δd .^ 2 .+ Δf .^ 2), label = "difference")
end
function plot_double_gaussian_fit_difference(plt, fid, run, βc; kws...)
    ups, f, Δf = histogram_jackknife_fit(fid, run, βc)
    ups, P, ΔP, covP, V, dS = probability_density(fid, run, βc)
    d, Δd = P * 6V, ΔP * 6V
    return plot!(plt, ups, d - f; ribbon = sqrt.(Δd .^ 2 .+ Δf .^ 2), label = "difference")
end
