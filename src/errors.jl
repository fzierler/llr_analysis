function mean_std_of_mean_cov(x0::AbstractArray; dims::Integer)
    N = size(x0, dims)
    x = dropdims(mean(x0; dims); dims)
    Δx = dropdims(std(x0; dims); dims) ./ sqrt(N)
    cov_x = cov(x0; dims) ./ sqrt(N)
    return x, Δx, cov_x
end
function mean_std_of_mean(x0::AbstractArray; dims::Integer)
    N = size(x0, dims)
    x = dropdims(mean(x0; dims); dims)
    Δx = dropdims(std(x0; dims); dims) ./ sqrt(N)
    return x, Δx
end
function mean_std_of_mean(x0::AbstractVector)
    N = length(x0)
    x = mean(x0)
    Δx = std(x0) ./ sqrt(N)
    return x, Δx
end
function apply_jackknife(obs::AbstractVector)
    N = length(obs)
    O = mean(obs)
    ΔO = sqrt(N - 1) * std(obs, corrected = false)
    return O, ΔO
end
function apply_jackknife(obs::AbstractArray; dims::Integer)
    N = size(obs)[dims]
    O = dropdims(mean(obs; dims); dims)
    ΔO = dropdims(sqrt(N - 1) * std(obs; dims, corrected = false); dims)
    return O, ΔO
end
