using Makie

abstract type Probe end

struct ProbeLong <: Probe
    num_channels::Int64
    spacing_column::Float64
    spacing_row::Float64
    length::Float64
    width::Float64
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, probe::ProbeLong)
    nrows = div(probe.num_channels,2)
    pos = Vector{Point2f0}(undef, probe.num_channels)
    for i in 1:probe.num_channels
        r = div(i-1,2)
        c = i-r*2-1
        pos[i] = Point2f0(c*103.0+6.0, r*20.0+6.0)
    end
    [PlotSpec(Scatter, pos;markersize=12, marker=:rect, markerspace=:data)]
end



