using LinearAlgebra

struct DirectionFiltered
    anglebins::AbstractVector{<:Real}
    weight::Dict{CartesianIndex{4}, Float64}
    occupancy::Dict{CartesianIndex{4}, Float64}
    index::Vector{CartesianIndex{5}}
end

function get_view_rate_map(dd::DirectionFiltered, nbins)
    X = zeros(nbins)
    Y = zeros(nbins)
    for k in keys(dd.weight)
        k1 = getindex(k,1)
        X[k1] += dd.weight[k]
        Y[k1] += dd.occupancy[k] 
    end
    X, Y
end

function get_view_rate_map(dd::DirectionFiltered, nbins, direction_ranges::UnitRange...)
    X = zeros(nbins, length(direction_ranges))
    Y = zeros(nbins, length(direction_ranges))
    for k in keys(dd.weight)
        k1 = getindex(k,1)
        k4 = getindex(k,4)
        for (jj,kr) in enumerate(direction_ranges)
            if k4 in kr
                X[k1,jj] += dd.weight[k]
                Y[k1,jj] += dd.occupancy[k] 
                break
            end
        end
    end
    X, Y
end

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
    # if the field was not exited at all, just use the last point
    if exited == 0
        exited = size(pos,2)
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

function DirectionFiltered(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, mm::SimpleMesh, idx::AbstractVector{<:Integer})
    nt = numtrials(qdata)
    θbins = range(-π, stop=π, length=24)
    gidx = CartesianIndex{5}[]
    weight = Dict{CartesianIndex{4}, Float64}()
    occupancy = Dict{CartesianIndex{4},Float64}()
    nn = zeros(Int64, nt)
    nm = zeros(Int64, nt)
    for i in 1:nt
        # find the first point at which the field is entered and when it is exited
        idx0,idx1 = get_direction(qdata.position[i][1:2,:], mm,idx)
        # debug: check the number of times the field is visited
        qmidx = jocc.index[i][vpvrp.placeviewidx[i]]
        nn[i] = length(filter(k->in(idx)(getindex(k,2)), qmidx))
        # TODO: Also get the gaze for these positions
        if idx1 >= idx0 > 0
            v = qdata.position[i][1:2,idx1] - qdata.position[i][1:2,idx0]
            θ = atan(v[2],v[1])
            l = searchsortedfirst(θbins, θ)
            vidx = findall(in(idx0:idx1), vpvrp.placeviewidx[i])
            # get the mesh bin indices
            qidx = jocc.index[i][vpvrp.placeviewidx[i][vidx]]
            # count the number of spikes
            # cc = length(vidx)
            # kk = filter(k->(k[2] in idx), qidx)
            for k in qidx
                ki = CartesianIndex(k[1], k[2], k[3], i)
                if ki in keys(jocc.weight)
                    nm[i] += 1
                    kk = CartesianIndex(k[1], k[2], k[3], l)
                    weight[kk] = get(weight, kk, 0.0)  + 1.0
                    occupancy[kk] = get(occupancy, kk, 0.0) + jocc.weight[ki]
                    push!(gidx, CartesianIndex(k[1], k[2], k[3], l, i))
                end
            end
        end
    end
    DirectionFiltered(θbins, weight, occupancy, gidx ), nn, nm
end

function get_directionality(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, rf::T) where T <: AbstractResponseFields
    mm = get_mesh(T, rf.args[:nrefinements])
    clusters = merge_fields(rf)
    nt = numtrials(qdata)
    λ = fill(NaN, nt, length(clusters))
    θ = fill(NaN, nt, length(clusters))
    gaze = Matrix{Matrix{Float64}}(undef, nt, length(clusters))
    for (ii,cluster) in enumerate(clusters)
        λ[:,ii], θ[:,ii],gg = get_directionality(qdata, vpvrp, mm, rf.binidx[cluster]) 
        for j in 1:nt
            if isassigned(gg, j, ii)
                gaze[j,ii] = gg[j,ii]
            end
        end
    end
    λ, θ, gaze
end

function get_direction_tuning(gidx::DirectionFiltered;smooth=false, α=0.1, niter=100)
    nc = length(gidx.index)
    X = zeros(length(gidx.anglebins),nc)
    Y = zeros(length(gidx.anglebins)nc)
    for j in 1:nc
        for (k,v) in gidx.occupancy[j]
            l = getindex(k,4)
            Y[l,j] += v
            if k in keys(gidx.weight[j])
                X[l,j] += gidx.weight[j][k]
            end
        end
    end
    if smooth
        Ls = get_normalize_laplacian(length(X))
        X = permutedims(laplace_smoothing(permutedims(X), Ls, α;niter=niter))
        Y = permutedims(laplace_smoothing(permutedims(Y), Ls, α;niter=niter))
    end
    X ./ Y
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

function plot_directional_view_field(dd::DirectionFiltered, args...;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_directional_view_field!(lg, dd, args...;kwargs...)
        link_cameras_lscene(fig)
        fig
    end
end

function plot_directional_view_field!(lg, dd::DirectionFiltered, place_field_idx::AbstractVector{<:Integer};floor_offset=-30)
    # south to north vs north to south
    # TODO: Should be tailored to each place field
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    X,Y = get_view_rate_map(dd, 1312,1:12, 13:24)
    Ls = get_normalize_laplacian(mm)
    Xs = laplace_smoothing(permutedims(X), Ls, 0.1;niter=100);
    Ys = laplace_smoothing(permutedims(Y), Ls, 0.1;niter=100);
    λ = Xs./Ys
    cr = extrema(filter(isfinite, λ))
    m_floor = Translate(0.0, 0.0, floor_offset)(Hippocampus.floor_topology3(;nrefinements=3))
    cm = Meshes.Point(mean(to.(centroid.(m_floor[place_field_idx])))...)
    cmp = Point3f(ustrip(cm.coords.x), ustrip(cm.coords.y), ustrip(cm.coords.z))
    Z = zeros(nelements(m_floor))
    Z[place_field_idx] .= 1.0
    Label(lg[1,1], "North → South", tellwidth=false)
    lscene1 = LScene(lg[2,1], show_axis=false)
    plotmesh!(lscene1, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene1, mm;color=Xs[1,:]./Ys[1,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene1, m_floor;color=Z)
    # indicate directionality
    arrows3d!(lscene1, cmp-Point3f(0.0,-2.5,0.0), Point3f([0.0, -5.0, 0.0]))

    Label(lg[1,2], "South → North",tellwidth=false)
    lscene2 = LScene(lg[2,2], show_axis=false)
    plotmesh!(lscene2, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene2, mm;color=Xs[2,:]./Ys[2,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene2, m_floor;color=Z)
    arrows3d!(lscene2, cmp-Point3f(0.0,2.5,0.0), Point3f([0.0, 5.0, 0.0]))

    # combined
    Xsc = laplace_smoothing(dropdims(sum(X,dims=2),dims=2), Ls, 0.1;niter=100);
    Ysc = laplace_smoothing(dropdims(sum(Y,dims=2),dims=2), Ls, 0.1;niter=100);
    λc = Xsc./Ysc
    Label(lg[1,3], "Combined", tellwidth=false)
    lscene3 = LScene(lg[2,3], show_axis=false)
    plotmesh!(lscene3, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene3, mm;color=λc, ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene3, m_floor;color=Z)

    # occupancy only 
    Label(lg[1,4], "Occupancy", tellwidth=false)
    lscene4 = LScene(lg[2,4], show_axis=false)
    plotmesh!(lscene4, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene4, mm;color=Ysc, ceiling_offset=10, floor_offset=-15)
    viz!(lscene4, m_floor;color=Z)


end