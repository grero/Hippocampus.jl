using LinearAlgebra

"""
    get_direction(pos::Matrix{<:Real}, mm::SimpleMesh, binidx::Vector{<:Integer})

Csompute the directional vector through the field characteried by `binidx`
"""
function get_direction(pos::Matrix{<:Real}, mm::SimpleMesh, binidx::AbstractVector{<:Integer})
    kn = KNearestSearch(mm,1)
    entered = 0 
    exited = 0
    for (ii,_pos) in enumerate(eachcol(pos))
        _idx,dd = searchdists(Meshes.Point(_pos...), kn)
        idx = first(_idx)
        if entered == 0
            if idx in binidx
                entered = ii 
            end
        end
        if entered > 0 && exited == 0
            if (!in(binidx))(idx)
                exited = ii
            end
        end
    end
    if entered > 0 && exited > 0
        v = pos[:,exited] - pos[:,entered]
        v ./= norm(v)
    end
    return entered,exited 
end

function get_direction(pos::Vector{Matrix}, args...)
    n = length(pos)
    p0 = zeros(2,n) 
    p1 = zeros(2,n)
    for (ii,_pos) in enumerate(pos)
        p0[:,ii], p1[:,ii] = get_direction(_pos, args...)
    end
    p0,p1
end

function get_direction(udata::UnityData, args...)
    nt = numtrials(udata) 
    p0 = zeros(2,nt)
    p1 = zeros(2,nt)
    for i in 1:nt
        tu,posx,posy,_ = get_trial(udata, i;trial_start=2)
        p0[:,i], p1[:,i] = get_direction(permutedims([posx posy]), args...)
    end
    p0,p1
end

function get_directionality(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, mm::SimpleMesh, idx::AbstractVector{<:Integer})
    nt = numtrials(qdata)
    λ = fill(NaN,nt)
    θ = fill(NaN,nt)
    for i in 1:nt
        idx0,idx1 = get_direction(qdata.position[i][1:2,:], mm,idx)
        if idx1 >= idx0 > 0
            v = qdata.position[i][1:2,idx1] - qdata.position[i][1:2,idx0]
            θ[i] = atan(v[2],v[1])
            vidx = findall(in(idx0:idx1), vpvrp.placeviewidx[i])
            # count the number of spikes
            cc = length(vidx)
            tt = qdata.timestamps[i][idx0:idx1]
            # we need to worry about gaps where
            dt = diff(tt)
            occ = sum(dt[dt.< 0.002])
            λ[i] = cc/occ
        end
    end
    λ,θ
end

function get_directionality(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
    clusters = merge_fields(rf)
    nt = numtrials(qdata)
    λ = fill(NaN, nt, length(clusters))
    θ = fill(NaN, nt, length(clusters))
    for (ii,cluster) in enumerate(clusters)
        λ[:,ii], θ[:,ii] = get_directionality(qdata, vpvrp, mm, rf.binidx[cluster]) 
    end
    λ, θ
end

function plot_field_directionality!(lg::GridLayout, λ::AbstractVector{<:Real}, θ::AbstractVector{<:Real};kwargs...)
    # compute circular mean
    qidx = isfinite.(λ)
    x = abs(sum(λ[qidx].*exp.(-im*θ[qidx]))/sum(λ[qidx]))
    θbins = range(-π, stop=π, length=24)
    h = fit(Histogram, θ, Weights(λ),θbins)
    ax = PolarAxis(lg[1,1],)
    hidedecorations!(ax)
    lines!(ax, θbins, [h.weights;h.weights[1]];linewidth=2.0, kwargs...)
    ax.thetaticklabelsvisible = false
    ax.rticksvisible = true
    ax.thetagridvisible=true
    ax.rgridvisible=true
    Label(lg[1,2], L"\mu_r=%$(round(x, sigdigits=2))", tellheight=false,rotation=-π/2)
end

function plot_field_directionality(λ::AbstractArray{<:Real}, θ::AbstractArray{<:Real}, args...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_field_directionality!(lg, λ, θ, args...)
        fig
    end
end

function plot_field_directionality!(lg, λ::Matrix{<:Real}, θ::Matrix{<:Real}, rf::T) where T <: AbstractResponseFields
    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    plot_response_fields!(lg1, rf)
    lg2 = GridLayout(lg[1,2])
    Label(lg2[1,1,TopLeft()], "B")
    lgi = [GridLayout(lg2[i,1]) for i in 1:size(θ,2)]
    colors = Makie.wong_colors()
    for (ii,_lg) in enumerate(lgi)
        plot_field_directionality!(_lg, λ[:,ii], θ[:,ii];color=colors[ii])
    end
    colsize!(lg, 1, Relative(0.7))
end