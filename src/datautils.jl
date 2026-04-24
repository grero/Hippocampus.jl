mutable struct RunCounter
    level::Int64
end

function (lc::RunCounter)(fname::String)
    lm = lc.level > 0
    lc.level -= 1
    return lm > 0
end

struct LevelCounter
    origin::String
    idx::Int64
    depth::Int64
end

function LevelCounter(origin::String, d)
    idx = findfirst(DPHT.levels.==origin)
    LevelCounter(origin, idx, d)
end

function (lc::LevelCounter)(fname::String,cwd=pwd())
    ll = DPHT.level(cwd)
    idx = findfirst(DPHT.levels.==ll)
    d = lc.idx-idx
    return d <= lc.depth 
end