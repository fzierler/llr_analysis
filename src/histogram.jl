# In python, we are using the following precision
# Mpmath settings:
#   mp.prec = 53                [default: 53]
#   mp.dps = 15                 [default: 15]
#   mp.trap_complex = False     [default: False]
setprecision(BigFloat, 53)

function log_partition_function(a, S, beta, ::Type{T} = BigFloat) where {T}
    # David uses a different sign for a
    dS = T(S[2]) - T(S[1])
    pi_exp = T(0)
    Z = T(0)
    exponent::Vector{T} = zeros(T, length(a))
    weight::Vector{T} = zeros(T, length(a))
    for (i, ai, Si) in zip(eachindex(a), a, S)
        # TODO: Use Kahan summation for a sums and cumulative sums
        A = T(beta) - T(ai)
        exponent[i] = T(pi_exp) + T(Si) * T(beta) - T(ai) * T(dS) / 2
        sinh_factor = iszero(A) ? dS / 2 : sinh(A * dS / 2) / A
        pi_exp -= ai * dS
        # TODO: perform weighted logsumexp here
        weight[i] = 2 * sinh_factor
    end
    Z = sum(exp.(exponent) .* weight)
    return log(Z)
end

"""
    energy_moment(S,a,β,N,[::Type{T}=Float64,::Type{U}=BigFloat])

    Calculate N-th cumulant of the energy from the d.o.s. coefficients `a`
    at the central energies `S` at the inverse coupling β.

    This code has been ported over from David Mason's implementation
    in python (10.5281/zenodo.13807993). It follows Eqs. (3.1.18) and
    (3.1.19) from David Masons's PhD thesis.

    By default, we use julia's BigFloat together with the double precision Float64
    datatype for floating point calculations. This can be overwritten by specifying
    the datatypes `T` and `U` to be used instead.

    Typically, it suffices to use `T=Float128` for all calculations in here as long
    as log(Z) is calculated in higher precision e.g. `U=BigFloat`.
"""
function energy_moment(S, a, β, N::Int, ::Type{T} = Float64, ::Type{U} = BigFloat, logZ = log_partition_function(a, S, β, U)) where {T, U}
    pi_exp::T = - T(logZ)
    full_exp::T = T(0)
    En::T = T(0)
    δS::T = T(S[2]) - T(S[1])
    for (Si, ai) in zip(S, a)
        A::T = - T(ai) + T(β)
        full_exp = 2 * exp(T(pi_exp) + T(β) * (T(Si) - T(δS) / 2) + T(A) * T(δS) / 2)
        for m in 0:(N)
            sinh_term::T = T(0)
            cosh_term::T = T(0)
            for j in 0:div(m, 2, RoundDown)
                sinh_term += (δS / 2)^(2j) * T(Si)^(m - 2j) / factorial(2j) / factorial(m - 2j)
            end
            if m > 0 # this condition skips terms proportional to 1/(-1)! → 1/Γ(0) → 0
                for j in 1:div(m, 2, RoundUp)
                    cosh_term += (δS / 2)^(2j - 1) * T(Si)^(m - 2j + 1) / factorial(2j - 1) / factorial(m - 2j + 1)
                end
            end
            sh = sinh(A * δS / 2)
            ch = cosh(A * δS / 2)
            Ap = factorial(N) * A^(m - N - 1) * (-1)^(N - m)
            En += full_exp * (sh * sinh_term + ch * cosh_term) * Ap
        end
        pi_exp -= ai * δS
    end
    return En
end

"""
    all_energy_moments(S,a,β,N,[::Type{T}=Float64,::Type{U}=BigFloat])

    Calculate all from the first to the N-th moment of the energy from the dos
    coefficients `a` at the central energies `S` at the inverse coupling β.

    Return as vector with `N` elements, with their indices corresponding
    to the moment of the plaquette
"""
function all_energy_moments(S, a, β, N::Int, ::Type{T} = Float64, ::Type{U} = BigFloat, logZ = log_partition_function(a, S, β, U)) where {T, U}
    pi_exp::T = - T(logZ)
    full_exp::T = T(0)
    En::Vector{T} = zeros(T, N)
    δS::T = T(S[2]) - T(S[1])
    for (Si, ai) in zip(S, a)
        A::T = - T(ai) + T(β)
        full_exp = 2 * exp(T(pi_exp) + T(β) * (T(Si) - T(δS) / 2) + T(A) * T(δS) / 2)
        for m in 0:(N)
            sinh_term::T = T(0)
            cosh_term::T = T(0)
            for j in 0:div(m, 2, RoundDown)
                sinh_term += (δS / 2)^(2j) * T(Si)^(m - 2j) / factorial(2j) / factorial(m - 2j)
            end
            if m > 0 # this condition skips terms proportional to 1/(-1)! → 1/Γ(0) → 0
                for j in 1:div(m, 2, RoundUp)
                    cosh_term += (δS / 2)^(2j - 1) * T(Si)^(m - 2j + 1) / factorial(2j - 1) / factorial(m - 2j + 1)
                end
            end
            sh = sinh(A * δS / 2)
            ch = cosh(A * δS / 2)
            factor = full_exp * (sh * sinh_term + ch * cosh_term)
            for i in max(1, m):N
                Ap = factorial(i) * A^(m - i - 1) * (-1)^(i - m)
                En[i] += factor * Ap
            end
        end
        pi_exp -= ai * δS
    end
    return En
end

_E_in_interval(E, S0, dS) = E >= (S0 - dS / 2) && E < (S0 + dS / 2)
function log_rho(E, S, dS, a; cumsum_a = cumsum(a))
    for i in eachindex(S, a)
        if _E_in_interval(E, S[i], dS)
            cs = i == 1 ? zero(dS) : cumsum_a[i - 1]
            log_ρ = a[i] * (S[i] - dS / 2 - E) - cs * dS
            return log_ρ
        end
    end
    return
end
function _set_up_histogram(fid, run)
    a, S_all = LLRParsing.a_vs_central_action_repeats(fid, run; ind = nothing)[1:2]
    Ns = read(fid[run], "Ns")::Int
    Nt = read(fid[run], "Nt")::Int
    S = unique(S_all)
    V = Nt * Ns^3
    return a::Matrix{Float64}, S::Vector{Float64}, Nt, Ns, V
end
function probability_density_repeats(fid, run, beta; kws...)
    a, S, Nt, Ns, V = _set_up_histogram(fid, run)
    return probability_density_repeats(a, S, beta, V; kws...)
end
function trapz(x, y)
    s = zero(eltype(x))
    @assert issorted(x)
    for i in eachindex(x)
        i == 1 && continue
        s += (y[i] + y[i - 1]) * (x[i] - x[i - 1]) / 2
    end
    return s
end
function probability_density_repeats(a, S, beta, V; nbins = length(S), normalize = false)
    dS = S[2] - S[1]
    E = range(minimum(S) + dS, maximum(S), length = nbins)
    ups = @. E / (6V)
    repeats = last(size(a)[2])
    P = zeros(nbins, repeats)

    for i in 1:repeats
        logZ = Float64(log_partition_function(a[:, i], S, beta))
        csa = cumsum(a[:, i])
        for j in 1:nbins
            log_ρ = log_rho(E[j], S, dS, a[:, i]; cumsum_a = csa)
            P[j, i] = exp(log_ρ + beta * E[j] - logZ)
        end
        if normalize
            #not really needed. ∫PdE is already approx. 1
            intP = trapz(E, P[:, i])
            @. P[:, i] = P[:, i] / intP
        end
    end
    return ups, P, V, dS
end
function probability_density(fid, run, beta; kws...)
    ups, prob, V, dS = probability_density_repeats(fid, run, beta; kws...)
    P, ΔP, covP = mean_std_of_mean_cov(prob, dims = 2)
    return ups, P, ΔP, covP, V, dS
end
function probability_density(a, S, beta, V; kws...)
    ups, prob, V, dS = probability_density_repeats(a, S, beta, V; kws...)
    P, ΔP, covP = mean_std_of_mean_cov(prob, dims = 2)
    return ups, P, ΔP, covP, V, dS
end
function plot_plaquette_histogram!(plt, fid, run, beta; kws...)
    a, S, Nt, Ns, V = _set_up_histogram(fid, run)
    ups, P, ΔP, covP, V, dS = probability_density(a, S, beta, V)
    label = LLRParsing.fancy_title(run)
    xlabel = L"u_p"
    ylabel = L"P_{\beta}(u_p)"
    plot!(plt; xlabel, ylabel, yticks = :none, left_margin = 5Plots.mm)
    plot!(plt, ups, P * 6V; label, ribbon = ΔP * 6V, kws...)
    # restric plotting range
    i1 = findfirst(x -> x > 1.0e-14, P * 6V)
    i2 = findlast(x -> x > 1.0e-10, P * 6V)
    plot!(plt, xlims = (ups[i1], ups[i2]))
    return plt
end
