using Makie
import Makie.SpecApi as S
using Clustering
using MultivariateStats
using LinearAlgebra
using LinearRegressionUtils

"""
A spatial representation of events
"""
struct SpatialRepresentation
    position::Vector{Vector{Point2f}}
    events::Vector{Vector{Float64}}
end

function SpatialRepresentation(spikes::Spiketrain, rp::RippleData, udata::UnityData)
    nt = numtrials(udata)
    position = Vector{Vector{Point2f}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    sp = spikes.timestamps/1000.0 #convert to seconds
    for i in 1:nt
        tp,posx,posy,_ = get_trial(udata,i)
        tp .-= tp[1]
        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])
        # align to trial start
        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        events[i] = sp_trial
        position[i] = Vector{Point2f}(undef, nspikes)
        for j in 1:nspikes
            k = searchsortedfirst(tp,sp_trial[j])
            if 0 < k <= length(posx)
                position[i][j] = Point2f(posx[k],posy[k])
            end
        end
    end
    SpatialRepresentation(position,events)
end

function SpatialRepresentation()
    sptrain = Spiketrain()
    rdata = cd(DPHT.process_level(level(RippleData))) do
        RippleData()
    end
    udata = cd(DPHT.process_level(level(UnityData))) do
        UnityData()
    end
    SpatialRepresentation(sptrain, rdata, udata)
end

numtrials(spr::SpatialRepresentation) = length(spr.position)

function Makie.convert_arguments(::Type{<:AbstractPlot}, spr::SpatialRepresentation)
    nt = numtrials(spr)
    points = Point2f[]
    for pos in spr.position
        append!(points, Point2f.(pos))
    end
    ax = S.Axis(plots=[S.Scatter(points)])
    S.GridLayout(ax)
end

function visualize!(lscene, spr::SpatialRepresentation;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0))
    nt = numtrials(spr)
    trial_events = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return spr.position[_trial.i]
        else
            return [Point2f(NaN)]
        end
    end
    scatter!(lscene, trial_events)
end

"""
Contains information about the total time spent in each spatial bin.
"""
struct SpatialOccupancy{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
end

function SpatialOccupancy(udata::UnityData, xbins::AbstractVector{T}, ybins::AbstractVector{T};trial_start=1) where T <: Real
    nt = numtrials(udata)
    weight = fill(0.0, length(xbins)-1, length(ybins)-1)
    for i in 1:nt
        tu, posx, posy, _ = get_trial(udata, i;trial_start=trial_start)
        for j in 2:length(tu)
            Δt = tu[j]-tu[j-1]
            xidx = searchsortedlast(xbins, posx[j-1])
            yidx = searchsortedlast(ybins, posy[j-1])
            if 0 < xidx <= size(weight,1) && 0 < yidx <= size(weight,2)
                weight[xidx,yidx] += Δt
            end
        end
    end
    SpatialOccupancy(xbins, ybins, weight)
end

function SpatialOccupancy(xbins,ybins=xbins;kwargs...)
    udata = UnityData()
    SpatialOccupancy(udata, xbins, ybins;kwargs...)
end

struct SpatialMap{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    weight::Matrix{T}
    occupancy::Matrix{T}
end

DPHT.level(::Type{SpatialMap}) = "cell"

function SpatialMap(spr::SpatialRepresentation, xbins::AbstractVector{T}, ybins::AbstractVector{T},spoc::SpatialOccupancy;kwargs...) where T <: Real
    spatial_count = fill(0.0, length(xbins)-1, length(ybins)-1)
    nt = numtrials(spr)
    for i in 1:nt
        position = spr.position[i]
        xpos = [pos[1] for pos in position]
        ypos = [pos[2] for pos in position]
        h = fit(Histogram, (xpos,ypos), (xbins, ybins))
        spatial_count .+= h.weights
    end
    SpatialMap(xbins,ybins, spatial_count, spoc.weight)
end

function SpatialMap(xbins, ybins;redo=false, do_save=false,kwargs...)
    sp = Spiketrain()
    rp = cd(DPHT.process_level(RippleData;kwargs...)) do
        RippleData()
    end
    udata = cd(DPHT.process_level(UnityData;kwargs...)) do
        UnityData()
    end
    spoc = SpatialOccupancy(udata, xbins, ybins)
    spr = SpatialRepresentation(sp, rp, udata;kwargs...);
    SpatialMap(spr, xbins, ybins, spoc)
end

function adaptive_smoothing(spm::SpatialMap, α=10000.0^2;filter_unoccupied=true)
    Z = adaptive_smoothing(spm.weight, spm.occupancy, α)
    if filter_unoccupied
        Z[spm.occupancy.==0] .= eltype(spm.weight)(NaN) 
    end
    Z
end

function compute_sic(spm::SpatialMap)

    x = spm.weight./spm.occupancy
    idx = isfinite.(x)
    p = spm.occupancy[idx]./sum(spm.occupancy[idx])
    r = sum(p.*x[idx])
    xr = x[idx]./r
    ll = log2.(xr)
    lidx = isfinite.(ll)
    sic = sum((p.*xr.*ll)[lidx])
end

function compute_entropy(sp::Union{SpatialMap, SpatialOccupancy})
    pp = sp.weight./sum(sp.weight)
    -sum(filter(isfinite, pp.*log2.(pp)))
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, spm::SpatialMap,args::NamedTuple=(;))
    default_arguments = Dict(:normalize=>true, :do_smooth=>true, :α=>10_000.0^2, :filter_gaze => true)
    for k in keys(args)
        v = args[k] 
        default_arguments[k] = v
    end
    X = copy(spm.weight)
    if default_arguments[:normalize]
        X = spm.weight./spm.occupancy
        label = "Firing rate [Hz]"
    else
        X = spm.weight
        label = "Spike count"
    end
    if default_arguments[:do_smooth]
        X = adaptive_smoothing(spm.weight, spm.occupancy, default_arguments[:α])
    end
    if default_arguments[:filter_gaze]
        X[spm.occupancy.==0] .= eltype(X)(NaN)
    end
    h = S.Heatmap(spm.xbins, spm.ybins, rotr90(X), colormap=:turbo)
    ax1 = S.Axis(plots=[h])
    ll = S.Colorbar(h,label=label)
    S.GridLayout([ax1 ll])
end

function Makie.convert_arguments(T::Type{<:AbstractPlot}, spm::SpatialMap,args::Vector{<:NamedTuple})
    a = [convert_arguments(T, spm, _arg) for _arg in args]
    S.GridLayout(a)
end

"""
Combine the spatial representions for each of the cells represented by `cellidirs` and use
those combined responses to regress to position in the maze
"""
function regress_space(celldirs::Vector{String};kwargs...)
    spr = map(celldirs) do celldir
       cd(celldir) do
            Hippocampus.SpatialRepresentation(;min_speed=0.3,trial_start=2)
        end
    end
    regress_space(spr;kwargs...)
end

function regress_space(spr::Vector{SpatialRepresentation{T1,T2}};n_spatial_clusters=256, kwargs...) where T1 <: Real where T2 <: Real
    X, position = get_population_representation(spr)
    km_results = kmeans(position, n_spatial_clusters)
    # sum up responses in each of the spatial bins returned by the kmean algorithm
    X2 = zeros(eltype(X), size(X,1), n_spatial_clusters)
    for (i,k) in enumerate(km_results.assignments)
       X2[:,k] .+= X[:,i]
    end
    
    # perform PCA to decorrelate the inputs
    pca = fit(PCA, X2)
    Z = predict(pca, X2)
    lq = LinearRegressionUtils.llsq_stats(permutedims(Z), permutedims(km_results.centers))
    lq, pca, km_results, X2,X, position
end

function plot_regression_reults(lq, km_results, position, X)
    # compute histogram of positions
    xbins = range(-12.5f0, stop=12.5f0, length=40)
    ybins = xbins
    h = fit(Histogram, (eachrow(position)...,), (xbins, ybins))
    
    # project X onto the regression space
    q,r = qr(lq.β[1:end-1,:])
    w = q[:,1:2]
    Y = w'*X
    # also do the actual regresson
    Yp = lq.β[1:end-1,:]'*X .+ lq.β[end,:]
    with_theme(plot_theme) do
        fig = Figure(size=(600,700))
        ax1 = Axis(fig[1,1])
        Label(fig[1,1,TopLeft()], "A")
        hh = heatmap!(ax1, xbins, ybins, h.weights)
        Colorbar(fig[1,2], hh, label="Count") 
        scatter!(ax1, Point2f.(eachcol(km_results.centers)), color=:white,
                markersize=5px)

        
        lg1 = GridLayout(fig[2,1:2])
        ax1_1 = Axis(lg1[1,1])
        Label(lg1[1,1,TopLeft()], "B")
        scatter!(ax1_1, km_results.centers[1,:], Y[1,:],color=km_results.centers[1,:],
                            colormap=:solar)
        ax1_1.xlabel = "True x-pos"
        ax1_1.ylabel = "Estimated x-pos"
        ax1_2 = Axis(lg1[1,2])
        Label(lg1[1,2,TopLeft()], "C")
        scatter!(ax1_2, km_results.centers[2,:], Y[2,:], color=km_results.centers[2,:],
                        colormap=:matter)
        ax1_2.xlabel = "True y-pos"
        ax1_2.ylabel = "Estimated y-pos"

        ax1_3 = Axis(lg1[1,3])
        Label(lg1[1,3,TopLeft()], "D")
        barplot!(ax1_3, [1], [lq.r²])
        ax1_3.ylabel = "r²"
        ax1_3.xticksvisible = false
        ax1_3.xticklabelsvisible = false
        colsize!(lg1, 3, Relative(0.15))

        lg2 = GridLayout(fig[3,1:2])
        ax2 = Axis(lg2[1,1])
        Label(lg2[1,1,TopLeft()],"E")
        sc = scatter!(ax2, Point2f.(eachcol(Y)), color=km_results.centers[1,:],colormap=:solar)
        #Colorbar(lg2[2,1], sc, label="x-pos", vertical=false, flipaxis=false)
        ax3 = Axis(lg2[1,2])
        Label(lg2[1,2,TopLeft()], "F")
        scy = scatter!(ax3, Point2f.(eachcol(Y)), color=km_results.centers[2,:], colormap=:matter)
        #Colorbar(lg2[2,2], scy, label="y-pos", vertical=false, flipaxis=false)
        #rowsize!(lg2, 2, Relative(0.25))
        fig
    end
end