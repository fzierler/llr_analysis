# For higher moments I use the didcated functions because they are faster
# I have checked that the numerical precision is sufficient
binder_cumulant(x) = 1 - moment(x, 4, 0.0) / mean(f -> f^2, x)^2 / 3
poly_susceptibility(x, Ns) = var(abs.(x)) * Ns^3
specific_heat_plaq(x, Ns, Nt) = var(x) * 6Ns^3 * Nt

# TODO: Give variance of the mean of those statistical uncertainties analytically
#       Bin data instead of separating the measurements by τint
function jackknife_resample_1d_reduction(x, f)
    resampled = similar(x)
    tmp = zeros(eltype(x), (length(x) - 1))
    for index in eachindex(x)
        for i in eachindex(x)
            index == i && continue
            j = i < index ? i : i - 1
            tmp[j] = x[i]
        end
        resampled[index] = f(tmp)
    end
    return apply_jackknife(resampled)
end
function std_observables(h5file, ens)
    f = h5open(h5file)
    plaq = f[ens]["plaquette"][]
    poly = f[ens]["polyakov_loop"][]
    Nt = f[ens]["Nt"][]
    Ns = f[ens]["Ns"][]
    return std_observables(plaq, poly, Nt, Ns)
end
function std_observables(plaq, poly, Nt, Ns)
    binder_plaq = binder_cumulant(plaq)
    sh_plaq = specific_heat_plaq(plaq, Ns, Nt)
    plaq_vev = mean(plaq)
    poly_sus = poly_susceptibility(poly, Ns)
    poly_vev = mean(abs, poly)

    # Estimate autocorrelation times:
    # To speed it up I first obtain the exponential autocorrelation time 'τexp' and then
    # compute the integrated autocorrelation time 'τint' using the Madras-Sokal windowing
    # technique on measurements seperated by 'τexp'.
    τexp_plaq = exponential_autocorrelation_time(plaq)
    step_plaq = Int(ceil(τexp_plaq))
    τint_plaq, Δτint_plaq = step_plaq .* madras_sokal_time(plaq[1:step_plaq:end])

    τexp_poly = exponential_autocorrelation_time(abs.(poly))
    step_poly = Int(ceil(τexp_poly / 4))
    τint_poly, Δτint_poly = step_poly .* madras_sokal_time(abs.(poly[1:step_poly:end]))

    # Determine uncertainties from data seperated by τint
    plaq_τ = plaq[1:Int(ceil(τint_plaq)):end]
    poly_τ = abs.(poly[1:Int(ceil(τint_poly)):end])

    Δbinder_plaq = jackknife_resample_1d_reduction(plaq_τ, x -> binder_cumulant(x))[2]
    Δsh_plaq =
        jackknife_resample_1d_reduction(plaq_τ, x -> specific_heat_plaq(x, Ns, Nt))[2]
    Δplaq_vev = std(plaq_τ) / sqrt(length(plaq_τ))
    Δpoly_sus = jackknife_resample_1d_reduction(poly_τ, x -> poly_susceptibility(x, Ns))[2]
    Δpoly_vev = std(poly_τ) / sqrt(length(poly_τ))

    # Obtain histograms like David did for further use
    hist = fit(Histogram, plaq, nbins = 100)
    hist = LinearAlgebra.normalize(hist, mode = :pdf)

    return binder_plaq,
        Δbinder_plaq,
        sh_plaq,
        Δsh_plaq,
        plaq_vev,
        Δplaq_vev,
        poly_sus,
        Δpoly_sus,
        poly_vev,
        Δpoly_vev
end
function write_all_std_observables_to_hdf5(h5file_in, h5file_out)
    fid = h5open(h5file_in)["ImportanceSampling"]
    return @showprogress for latt in keys(fid)
        for β in keys(fid[latt])
            plaq_vec = fid[joinpath(latt, β, "plaquette")][]
            poly_vec = fid[joinpath(latt, β, "polyakov_loop")][]
            Nt = fid[joinpath(latt, β, "Nt")][]
            Ns = fid[joinpath(latt, β, "Ns")][]
            beta = fid[joinpath(latt, β, "beta")][]

            bc, Δbc, sh, Δsh, plaq, Δplaq, psus, Δpsus, poly, Δpoly =
                std_observables(plaq_vec, poly_vec, Nt, Ns)

            # Save newly computed data to file
            h5write(h5file_out, joinpath("ImportanceSampling", latt, β, "plaquette"), plaq)
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "Delta_plaquette"),
                Δplaq,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "polyakov_loop"),
                poly,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "Delta_polyakov_loop"),
                Δpoly,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "binder_cumulant"),
                bc,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "poly_susceptibility"),
                psus,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "specific_heat_plaq"),
                sh,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "Delta_binder_cumulant"),
                Δbc,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "Delta_poly_susceptibility"),
                Δpsus,
            )
            h5write(
                h5file_out,
                joinpath("ImportanceSampling", latt, β, "Delta_specific_heat_plaq"),
                Δsh,
            )
            # Also save basic lattice parameters
            h5write(h5file_out, joinpath("ImportanceSampling", latt, β, "Nt"), Nt)
            h5write(h5file_out, joinpath("ImportanceSampling", latt, β, "Ns"), Ns)
            h5write(h5file_out, joinpath("ImportanceSampling", latt, β, "beta"), beta)

        end
    end
end
