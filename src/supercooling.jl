function supercooling(h5dset, run)
    a, S, Nt, Ns, V = LLRParsing._set_up_histogram(h5dset, run)
    a_jk = jackknife_resamples(a)

    t1_jk = zeros(size(a_jk)[2])
    t2_jk = zeros(size(a_jk)[2])
    tc_jk = zeros(size(a_jk)[2])

    for i in axes(a_jk, 2)
        ai = a_jk[:, i:i]
        βc = LLRParsing.beta_at_equal_heights(ai, S, V)
        pks = only(findmaxima(vec(ai), 5).indices)
        mns = only(findminima(vec(ai), 5).indices)
        t1_jk[i] = 1 / ai[pks]
        t2_jk[i] = 1 / ai[mns]
        tc_jk[i] = 1 / βc
    end

    t1, Δt1 = apply_jackknife(t1_jk)
    t2, Δt2 = apply_jackknife(t2_jk)
    tc, Δtc = apply_jackknife(tc_jk)
    return t1, Δt1, t2, Δt2, tc, Δtc
end
