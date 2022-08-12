function zscore_smooth_idxs(Y, lag, threshold, influence)
    N = length(Y)
    smooth = ones(Bool, N)
    filteredY = copy(Y)
    avgFilter = zeros(N)
    stdFilter = zeros(N)
    avgFilter[lag - 1] = mean(Y[1:lag])
    stdFilter[lag - 1] = std(Y[1:lag])
    @inbounds for i in range(lag, stop=N-1)
        if abs(Y[i] - avgFilter[i - 1]) > threshold  * stdFilter[i - 1]
            smooth[i] = false
            filteredY[i] = influence*Y[i] + (1 - influence)*filteredY[i - 1]
        else
            filteredY[i] = Y[i]
        end
        avgFilter[i] = mean(filteredY[i-lag+1:i])
        stdFilter[i] = std(filteredY[i-lag+1:i])
    end

    return smooth
end

function zscore_idxs(Y, threshold)
    N = length(Y)
    σ = std(Y)
    μ = mean(Y)
    lb = μ - σ * threshold
    ub = μ + σ * threshold
    (x -> lb < x < ub).(Y)
end

function frechet(A, B)
    C = pairwise(Euclidean(eps(eltype(A))), A, B, dims=2)
    for i in 2:size(A)[2]
        C[i, 1] = max(C[i-1, 1], C[i, 1])
    end

    for j in 2:size(B)[2]
        C[1, j] = max(C[1, j-1], C[1, j])
        for i in 2:size(A)[2]
            C[i, j] = max(min(C[i-1, j-1], C[i, j-1], C[i-1, j]), C[i, j])
        end
    end

    C[end, end]
end

function lsqfitting(X::Vector{T}, Y::Vector{T}) where {T<:Real}
    y0, y1 = minimum(Y), maximum(Y)
    x0, x1 = minimum(X), maximum(X)

    @variables x[1:length(X)] k n
    x = collect(x)
    hillexpr = y0 .+ (y1 - y0) .* k^n ./ (x.^n .+ k^n)
    j = Symbolics.jacobian(hillexpr, [k, n])
    J, _ = build_function(j, x, [k, n], expression=Val{false})

    @. F(x, p) = y0 + (y1 - y0) * p[1]^p[2] / (x^p[2] + p[1]^p[2])
    P = [(x1 + x0) / 2, one(T)]
    opts = Dict(
        :upper => [T(Inf), T(Inf)],
        :lower => [zero(T), zero(T)],
        :show_trace => false
    )
    k, n = curve_fit(F, J, X, Y, P; opts...).param
    y0, y1, k, n
end

function context_groupings(df; backbone=missing, strain=missing)
    df = copy(df)

    if !(strain === missing)
        filter!(r -> r.sA == r.sB == strain, df)
    end

    if !(backbone === missing)
        filter!(r -> r.bA == r.bB == backbone, df)
    end
    
    groups = groupby(df, [:pA, :pB])
    sort(combine(groups, :c => any, renamecols=false), [:pA, :pB])
end

function compatibilitymatrix(df)
    x = unique(df.pA)
    y = unique(df.pB)
    z = Matrix{String63}(undef, length(x), length(y))

    for i in 1:length(x)
        for j in 1:length(y)
            if df.c[(df.pA .== x[i]) .& (df.pB .== y[j])][1]
                z[i, j] = "Compatible"
            else
                z[i, j] = "Incompatible"
            end
        end
    end
    x, y, z
end
