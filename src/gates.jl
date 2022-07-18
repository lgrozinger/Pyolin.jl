function keep!(S::FlowSample, idxs)
    for channel in keys(S)
        S.data[channel] = S[channel][idxs]
    end
    S.params["\$TOT"] = string(length(idxs))
    S
end

function apply!(S::FlowSample, f, channels...)
    for channel in channels
        S.data[channel] = S.(S[channel])
    end
    S
end

Window!(start::Int, finish::Int) = S -> keep!(S, start:finish)
Head!(n::Int) = Window!(1, n)
Tail!(n::Int) = S -> Window!(length(S[CHANNEL]) - n)(S)
Filter!(f, channels...) = S -> keep!(S, f.([S[ch] for ch in channels]...))
NonZero!(channels...) = Filter!((xs...) -> all(xs .> 0), channels...)
Linear!(m, c, x, y) = Filter!((a, b) -> b < m * a + c, x, y)
Apply!(f, channels...) = S -> apply!(S, f, channels...)
Shift!(x, channels...) = Apply!(y -> x + y, channels...)
Scale!(x, channels...) = Apply!(y -> x * y, channels...)

function Spike!(channel, lag, threshold, influence)
    function (S)
        idxs = zscore_smooth_idxs(S[channel], lag, threshold, influence)
        keep!(S, idxs)
    end
end

function Zscore!(channel, threshold)
    S -> keep!(S, zscore_idxs(S[channel], threshold))
end

ZeroMin!(channel) = S -> Shift!(minimum(S[channel]), channel)

function Total!(S)
    tot = S["B1-A"] .+ S["B2-A"] .+ S["V1-A"] .+ S["V2-A"] .+ S["Y1-A"] .+ S["Y2-A"] .+ S["Y3-A"] .+ S["Y4-A"]
    scat = S["FSC-A"] .+ S["SSC-A"]
    keep!(S, tot .< scat)
end

function RPUConvert!(channel, auto, standard)
    yfp₀ = auto
    yfpₛ = standard
    function (S)
        yfp = median(S[CHANNEL])
        Scale!((yfp - yfp₀) / (yfp * (yfpₛ - yfp₀)), channel)
    end
end
