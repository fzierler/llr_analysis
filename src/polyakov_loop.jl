function is_fixed_a_measured(h5, ens)
    repeats = read(h5[ens],"repeats")
    # Read all central energies and coefficients a_n
    # If no fixed a-calculation was done, then an empty array will be returned 
    S0 = [ read(h5["$ens/$r"],"S0_fxa") for r in repeats]
    an = [ read(h5["$ens/$r"],"an_fxa") for r in repeats]
    # detrmines if we have done the measurement everywhere
    if any(isempty,S0) || any(isempty,an) 
        return false
    end
    return true
end
function read_fixed_data(h5,ens)
    repeats = read(h5[ens],"repeats")
    Nrep = read(h5[ens],"N_repeats")
    Nint = read(h5[ens],"N_replicas")
    # If we proceed with the plot then we have performed the measurements 
    # everywhere. Thus all entries of  S0 are idential and we can pick
    # one representative.
    # For the an's we have different values for every repeat and we rehape
    # the data into an array of size (N_intervals,N_repeats)
    S0 = [ read(h5["$ens/$r"],"S0_fxa") for r in repeats]
    an = [ read(h5["$ens/$r"],"an_fxa") for r in repeats]
    S0 = first(S0)
    an = hcat(an...)

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
    return an, S0, E, poly
end
# S0 is the same as Ek in David's code
# an is the same as -a in David's code
# E  is the sames as S in David's code 
# TODO: Look at performance eventually
function logZ_fixed_a(E, S0, an, β)
    dS = S0[2] - S0[1]
    # Determine the largest possible exponent for the first energy interval
    # (It will only be used to improve numerical stability)
    # This expression matches David's code
    lenE = length(E[1,:,:])
    log_ρ = LLRParsing.log_rho(S0[1], S0, dS, an)
    VEV_exp_0 = @. ( -an[1] + β) * E[1,:,:] + an[1]*S0[1] + log_ρ - log(lenE) + log(dS)
    vmax = maximum(VEV_exp_0)
    # Now add ap all contributions from every energy interval
    Z = zeros(length(S0))
    for i in eachindex(S0)
        log_ρ = LLRParsing.log_rho(S0[i], S0, dS, an)
        VEV_exp = @. (-an[i] + β) * E[i,:,:] + an[i]*S0[i] + log_ρ - log(lenE) + log(dS) - vmax
        Z[i] = sum(exp,VEV_exp)
    end
    logZ = vmax + log(sum(Z))
    # everything matches up to here
    return logZ
end
function polyakov_loop_fixed_a(E, S0, an, β, poly; f=abs)
    dS = S0[2] - S0[1]
    obs = zeros(length(S0))
    log_Z = logZ_fixed_a(E, S0, an, β)
    tmp = similar(E[1,:,:])

    for i in eachindex(S0)
        log_ρ = LLRParsing.log_rho(S0[i], S0, dS, an)
        @. tmp = (β * E[i,:,:]) - an[i]*(E[i,:,:] - S0[i]) + log_ρ - log_Z
        @. tmp = dS * f(poly[i,:,:]) * exp(tmp)
        obs[i] = mean(tmp)
    end
    res = sum(obs)
    return res
end
