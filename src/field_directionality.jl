using LinearAlgebra

struct DirectionFiltered
    anglebins::AbstractVector{<:Real}
    weight::Vector{Dict{CartesianIndex{4}, Float64}}
    occupancy::Vector{Dict{CartesianIndex{4}, Float64}}
    index::Vector{Vector{CartesianIndex{5}}}
end

function issignificant_old(gidx::DirectionFiltered;smooth=true, α=0.1, niter=1,pv_threshold=0.01)
    λ = get_direction_tuning(gidx;smooth=smooth, niter=niter, α=α) 
    z = λ.*exp.(gidx.anglebins*im)
    pv = zeros(size(z,2))
    for i in axes(z,2)
        pv[i] = pvalue(RayleighTest(filter(x->abs(x)>0, z[:,i])))
    end
    pv .< pv_threshold
end

function get_pvalue(gidx::DirectionFiltered;smooth=true, α=0.1, niter=1,nruns=1000, kwargs...)
    μr0,ϕ0 = get_directional_tuning_strength(gidx;smooth=smooth, α=α,niter=niter,kwargs...)
    μr = zeros(length(μr0),nruns)
    for i in 1:nruns
        μr[:,i],_ = get_directional_tuning_strength(gidx;do_shuffle=true, smooth=smooth, α=α, niter=niter, kwargs...)
    end
    pv = zeros(length(μr0))
    for i in 1:length(pv)
        if sum(isfinite.(μr[i,:])) > 20
            k = Makie.KernelDensity.kde(μr[i,:],boundary=(0,1))
            xx = sort(μr[i,:])
            cc = sum(pdf(k,xx[xx.<=μr0[i]]))/sum(pdf(k,xx))
            pv[i] = 1-cc
        end
    end
    pv
end

function issignificant(gidx::DirectionFiltered;pv_threshold=0.01,kwargs...)
    pv = get_pvalue(gidx;kwargs...)
    pv .< pv_threshold
end

function get_directional_tuning_strength(gidx;idx=1:length(gidx.anglebins), kwargs...)
    λ = get_direction_tuning(gidx;idx=idx, kwargs...)
    get_directional_tuning_strength(λ, gidx.anglebins)
end

function get_directional_tuning_strength(λ::Matrix{<:Real}, anglebins::AbstractVector{<:Real})
    z = λ.*exp.(anglebins*im)
    μr = zeros(size(λ,2))
    ϕ = zeros(size(λ,2))
    for i in axes(λ,2)
        _z = filter(x->abs(x)>0, z[:,i])
        zm = sum(_z)/sum(abs.(_z))
        ϕ[i] = angle(zm)
        μr[i] = abs(zm)
    end
    μr,ϕ
end

abstract type AbstractDirection end
struct East <: AbstractDirection end
Base.string(::Type{East}) = "East"
struct West <: AbstractDirection end
Base.string(::Type{West}) = "West"
struct North <: AbstractDirection end
Base.string(::Type{North}) = "North"
struct South <: AbstractDirection end
Base.string(::Type{South}) = "South"

function get_direction(gidx::DirectionFiltered, ::Type{East};Δθ=π/3)
    union(findall(gidx.anglebins .> π - Δθ/2), findall(gidx.anglebins .< -π+Δθ/2))
end

function get_arrow(mm::SimpleMesh, ::Type{East})
    (-2.5, 0.0,0.0), (5.0, 0.0, 0.0)
end

function get_direction(gidx, ::Type{West};Δθ=π/3)
    union(findall(0 .<= gidx.anglebins .< Δθ/2), findall(-Δθ/2 .<= gidx.anglebins .< 0))
end

function get_arrow(mm::SimpleMesh, ::Type{West})
    (2.5, 0.0,0.0), (-5.0, 0.0, 0.0)
end

function get_direction(gidx, ::Type{North};Δθ=π/3)
    findall( π/2 - Δθ/2 .<= gidx.anglebins .< π/2 + Δθ/2)
end

function get_arrow(mm::SimpleMesh, ::Type{North})
    (0.0, -2.5, 0.0), (0.0, 5.0, 0.0)
end

function get_direction(gidx, ::Type{South};Δθ=π/3)
    findall( -π/2 - Δθ/2 .< gidx.anglebins .<= -π/2 + Δθ/2)
end

function get_arrow(mm::SimpleMesh, ::Type{South})
    (0.0, 2.5, 0.0), (0.0, -5.0, 0.0)
end



function get_view_rate_map(dd::DirectionFiltered, nbins)
    X = zeros(nbins)
    Y = zeros(nbins)
    for k in keys(dd.occupancy)
        k1 = getindex(k,1)
        Y[k1] += dd.occupancy[k] 
        X[k1] += get(dd.weight, k, 0.0)
    end
    X, Y
end

function get_view_rate_map(dd::DirectionFiltered, ii::Integer, nbins, directions::Vector{Vector{Int64}})
    nq = length(directions)
    X = zeros(nbins, nq)
    Y = zeros(nbins, nq)
    for k in keys(dd.occupancy[ii])
        k1 = getindex(k,1)
        k4 = getindex(k,4)
        for (jj,kr) in enumerate(directions)
            if k4 in kr
                Y[k1,jj] += dd.occupancy[ii][k] 
                X[k1,jj] += get(dd.weight[ii],k,0.0)
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
        exited = size(pos,2)-1
    end
    if entered > 0 && exited > 0
        v = pos[:,exited] - pos[:,entered]
        v ./= norm(v)
    end
    return entered,exited 
end

function get_direction(pos::Vector{<:Matrix{<:Real}}, args...)
    n = length(pos)
    idx0 = fill(0,n) 
    idx1 = fill(0,n) 
    for (ii,_pos) in enumerate(pos)
        idx0[ii], idx1[ii] = get_direction(_pos[1:2,:], args...)
    end
    idx0,idx1
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

function DirectionFiltered(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, mm::SimpleMesh, rf::SpatialResponseFields;trial_start=1,kwargs...)
    nt = numtrials(qdata)
    θbins = range(-π, stop=π, length=24)
    nn = zeros(Int64, nt)
    nm = zeros(Int64, nt)
    clusters = merge_fields(mm, rf.binidx)
    nclusters = Hippocampus.get_num_fields(rf)
    cidx = findall(dropdims(mean(nclusters,dims=2),dims=2) .< 0.001)
    gidx = Vector{Vector{CartesianIndex{5}}}(undef, length(cidx))
    weight = Vector{Dict{CartesianIndex{4}, Float64}}(undef, length(cidx))
    occupancy = Vector{Dict{CartesianIndex{4},Float64}}(undef, length(cidx))
    for i in 1:length(cidx)
        gidx[i] = CartesianIndex{5}[]
        weight[i] = Dict{CartesianIndex{4}, Float64}()
        occupancy[i] = Dict{CartesianIndex{4}, Float64}()
    end
    for i in 1:nt
        for (ll,idx) in enumerate(clusters[cidx])
            # find the first point at which the field is entered and when it is exited
            tg,gaze,pos, fixmask,fo = get_trial(qdata,i;trial_start=trial_start);
            idx0,idx1 = get_direction(pos[1:2,:], mm,rf.binidx[idx])
            # debug: check the number of times the field is visited
            qmidx = jocc.index[i][vpvrp.placeviewidx[i]]
            nn[i] = length(filter(k->in(idx)(getindex(k,2)), qmidx))
            # TODO: Also get the gaze for these positions
        
            if idx1 >= idx0 > 0
                v = qdata.position[i][1:2,idx1] - qdata.position[i][1:2,idx0]
                θ = atan(v[2],v[1])
                l = searchsortedfirst(θbins, θ)
                # get the mesh bin indices
                # we actually do not want this; we want all the bins, not just the ones with spikes
                #this is the subset with spikes

                # occupancy
                # @show length(jocc.index[i]), idx0, idx1
                qpidx = jocc.index[i][idx0:idx1]
                tt = qdata.timestamps[i]
                for (jj,qp) in zip(idx0:idx1, qpidx)
                    if qp == CartesianIndex(0,0,0)
                        # skip invalid bins
                        continue
                    end
                    Δt = tt[jj+1]-tt[jj]
                    kk = CartesianIndex(qp[1], qp[2], qp[3], l)
                    occupancy[ll][kk] = get(occupancy[ll], kk, 0.0) + Δt
                    vv = findfirst(vpvrp.placeviewidx[i].==jj)
                    if vv !== nothing
                        weight[ll][kk] = get(weight[ll], kk, 0.0) + 1.0
                    end
                    push!(gidx[ll], CartesianIndex(qp[1], qp[2], qp[3], l, i))
                end
            end
        end
    end
    DirectionFiltered(θbins, weight, occupancy, gidx )
end

function DirectionFiltered(;redo=fname->false, do_save=true,kwargs...)
    h = process_kwargs(UnityRaytraceData;kwargs...)
    h = process_kwargs(ViewAndPlaceRepresentationNew,h;kwargs...)
    h = process_kwargs(JointOccupancy,h;kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    fname = "direction_filtered_placefields.jld2"
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(DirectionFiltered,fname)
    else
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        sessiondir = DPHT.get_level_path("session")
        qdata, jocc = cd(sessiondir) do
            qdata = UnityRaytraceData(;kwargs...)
            jocc = JointOccupancy(;kwargs...)
            # TODO: Also deal with other filtering here
            qdata, jocc
        end
        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields, 10_000;kwargs...)
        obj = DirectionFiltered(qdata, vpvrp, jocc, m_floor, rf_spatial;kwargs...)
        if do_save
            save_jld2(obj, fname)
        end
    end
    obj
end

function get_directionality(qdata::UnityRaytraceData, vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, rf::T) where T <: AbstractResponseFields
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

function get_direction_tuning(gidx::DirectionFiltered;idx=1:length(gidx.anglebins), do_shuffle=false, smooth=false, α=0.1, niter=1, spikes_only=false, occupancy_only=false)
    nc = length(gidx.index)
    X = zeros(length(gidx.anglebins),nc)
    Y = zeros(length(gidx.anglebins),nc)
    for j in 1:nc
        for (k,v) in gidx.occupancy[j]
            if do_shuffle
                l = rand(1:length(gidx.anglebins))
            else
                l = getindex(k,4)
            end
            Y[l,j] += v
            if k in keys(gidx.weight[j])
                X[l,j] += gidx.weight[j][k]
            end
        end
    end
    if smooth
        Ls = get_normalize_laplacian(size(X,1))
        X = permutedims(laplace_smoothing(permutedims(X[idx,:]), Ls, α;niter=niter))
        Y = permutedims(laplace_smoothing(permutedims(Y[idx,:]), Ls, α;niter=niter))
    else
        X = X[idx,:]
        Y = Y[idx,:]
    end
    if spikes_only
        return X
    elseif occupancy_only
        return Y
    end
    X ./ Y
end

"""
Return the bins in `mm` that wall within a 60 degree wedge centered on `pos`
"""
function get_view_bins(mm::SimpleMesh, pos, ϕ;Δϕ=π/3)
    midx = Int64[]
    # elevation in 120 degree wedge centered on π/2
    for ξ in range(π/2-π/3, stop=π/2+π/3, length=15)
        ϕ0 = ϕ-Δϕ/2
        ϕ1 = ϕ+Δϕ/2
        while ϕ0 < ϕ1 
            r0 = Ray(Tuple(pos), (sin(ξ)*cos(ϕ0),sin(ξ)*sin(ϕ0), cos(ξ)))
            # find all points where the ray intersects the maze
            iq0 = [Meshes.intersect(r0, m) for m in mm] 
            qidx0 = findall(iq0.!==nothing)
            idx0 = qidx0[argmin(norm.(centroid.(mm[qidx0]) .- pos))] 
            if !(idx0 in midx)
                push!(midx, idx0)
            end
            ϕ0 += Δϕ/15 # 10 steps
        end
    end
    midx
end

function get_view_bins(anglebins::AbstractVector{<:Real}, mm::SimpleMesh, m_floor::SimpleMesh, place_dir_idx::Vector{Tuple{Int64, Int64}};z=0.5, kwargs...)
    left_gaze_place_dir_idx = NTuple{3, Int64}[]
    right_gaze_place_dir_idx = NTuple{3, Int64}[]
    for (pidx, lidx) in place_dir_idx
        # player position
        pos = centroid.(m_floor[pidx]) + Meshes.Vec(0.0, 0.0, z)
        ϕ = anglebins[lidx]

        midx_left = get_view_bins(mm, pos, ϕ+π/2)
        for m in unique(midx_left)
            cc = (m,pidx,lidx)
            if !(cc in left_gaze_place_dir_idx)
                push!(left_gaze_place_dir_idx, (m, pidx, lidx))
            end
        end
        midx_right = get_view_bins(mm, pos, ϕ-π/2)
        for m in unique(midx_right)
            cc = (m,pidx,lidx)
            if !(cc in right_gaze_place_dir_idx)
                push!(right_gaze_place_dir_idx, (m, pidx, lidx))
            end
        end
    end 
    left_gaze_place_dir_idx, right_gaze_place_dir_idx
end

function get_egocentric_view_tuning(gidx, mm::SimpleMesh, m_floor::SimpleMesh;kwargs...)
    λright = Vector{Vector{Float64}}(undef, length(gidx.occupancy))
    λleft = Vector{Vector{Float64}}(undef, length(gidx.occupancy))
    for i in 1:length(gidx.occupancy)
        # get the place bins and traversal directions for this field
        vq = [(k[2], k[4]) for k in keys(gidx.occupancy[i])]
        left_gaze_place_dir_idx, right_gaze_place_dir_idx = get_view_bins(gidx.anglebins, mm, m_floor, unique(vq))
        # seperate the occupancy and weight into left and right
        λleft[i] = Float64[]
        λright[i] = Float64[]
        # create a view map first, by conditioning on place and direction
        # then sum up bins according to left_.. and right_.. above
        for (kk,vv) in gidx.occupancy[i]
            if (kk[1], kk[2], kk[4]) in left_gaze_place_dir_idx
                # TODO: Allow smoothing here
                #       This is the view map, 
                λ = get(gidx.weight[i], kk, 0.0)/vv
                push!(λleft[i], λ)
            elseif (kk[1], kk[2], kk[4]) in right_gaze_place_dir_idx
                λ = get(gidx.weight[i], kk, 0.0)/vv
                push!(λright[i], λ)
            end
        end
    end
    λleft, λright
end

"""
Estimate the preference for left vs right gaze relative to the direction of traverersal through the field
"""
function get_egocentric_gaze(gidx;pv_threshold=0.01,z=0.5)
    m_floor = (Hippocampus.floor_topology3(;nrefinements=3))
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    # this only makes sense for a direction field, so test that first
    # actually, we can still get egocentric coding even without directionality
    # the cell could be coding left/right gaze for any direction of traversal
    μr,ϕ = get_directional_tuning_strength(gidx;do_shuffle=false, smooth=true,niter=1)
    midx = Vector{Vector{CartesianIndex{6}}}(undef, length(μr))
    for i in 1:length(μr)
        midx[i] = CartesianIndex{6}[]
        for (j,(θ ,direction)) in enumerate(zip([π, π/2, 0.0, -π/2], [West, North, East, South]))
                # traversal direction
            aidx = get_direction(gidx, direction)
            for k in gidx.index[i]
                if k[4] in aidx
                    # the player is located z distance from the floor
                    pos = centroid(m_floor[k[2]]) + Meshes.Vec(0.0, 0.0, z)
                    # in the correct direction
                    for (kk,Δϕ) in enumerate([π/2, -π/2])
                        # TODO: Do we need to use multiple rays here; 
                        ϕ0 = θ+Δϕ-π/6
                        ϕ1 = θ+Δϕ+π/6
                        # TODO Also go vertical
                        for ξ in range(-π/6, stop=π/6, length=10)
                            while ϕ0 < ϕ1 
                                r0 = Ray(Tuple(pos), (sin(ξ)*cos(ϕ0),sin(ξ)*sin(ϕ0), cos(ξ)))
                                # find all points where the ray intersects the maze
                                iq0 = [Meshes.intersect(r0, m) for m in mm] 
                                qidx0 = findall(iq0.!==nothing)
                                idx0 = qidx0[argmin(norm.(centroid.(mm[qidx0]) .- pos))] 
                                cc = CartesianIndex(idx0, k[2], k[3], k[4], kk, k[5] )
                                if !(cc in mix[i])
                                    push!(midx[i], cc)
                                end
                                ϕ0 += π/30 # 10 steps
                            end
                        end
                    end

                end
                
            end
        end
    end
    midx
end


## plots

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

function plot_field_directionality(args...)
    with_theme(plot_theme) do
        fig = Figure(size=(700,400))
        lg = GridLayout(fig[1,1])
        plot_field_directionality!(lg,args...)
        fig
    end
end

function plot_field_directionality!(lg, gidx::DirectionFiltered, rf::T) where T <: AbstractResponseFields
    λ = get_direction_tuning(gidx)
    plot_field_directionality!(lg, λ, gidx.anglebins, rf)
end

function plot_field_directionality!(lg, λ::Matrix{<:Real}, θ::AbstractVector{<:Real}, rf::T) where T <: AbstractResponseFields
    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    plot_response_fields!(lg1, rf)
    lg2 = GridLayout(lg[1,2])
    Label(lg2[1,1,TopLeft()], "B")
    lgi = [GridLayout(lg2[i,1]) for i in 1:size(θ,2)]
    colors = Makie.wong_colors()
    for (ii,_lg) in enumerate(lgi)
        plot_field_directionality!(_lg, λ[:,ii], θ;color=colors[ii])
    end
    colsize!(lg, 1, Relative(0.7))
end

function plot_field_directionality!(lg, gidx::DirectionFiltered, rf::SpatialResponseFields;kwargs...)
    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    plot_response_fields!(lg1, rf)
    lg2 = GridLayout(lg[1,2])
    plot_directional_tuning!(lg2, gidx;kwargs...)
    Label(lg2[0,1,TopLeft()], "B")
    Label(lg2[0,1, Top()], "Firing rate", tellwidth=false, fontsize=14)
    rowsize!(lg2, 0, 5)

    lg3 = GridLayout(lg[1,3])
    plot_directional_tuning!(lg3, gidx;occupancy_only=true, kwargs...)
    Label(lg3[0,1,TopLeft()], "C")
    Label(lg3[0,1, Top()], "Occupancy", tellwidth=false, fontsize=14)
    rowsize!(lg3, 0, 5)
    colsize!(lg, 1, Relative(0.6))
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

function plot_directional_view_field!(lg, dd::DirectionFiltered, place_field_idx::AbstractVector{<:Integer},idx::Integer;floor_offset=-30, directions=[North, South])
    # south to north vs north to south
    # TODO: Should be tailored to each place field
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    binranges = [get_direction(dd, d) for d in directions]
    X,Y = get_view_rate_map(dd, idx, 1312,binranges)
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
   # Label(lg[1,1], "$(directions[1]) → $(directions[2])", tellwidth=false)
    lscene1 = LScene(lg[2,1], show_axis=false)
    plotmesh!(lscene1, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene1, mm;color=Xs[1,:]./Ys[1,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene1, m_floor;color=Z)
    # indicate directionality
    a,b = get_arrow(mm, directions[2])
    arrows3d!(lscene1, cmp-Point3f(a), Point3f(b))

   # Label(lg[1,2], "$(directions[2]) → $(directions[1])",tellwidth=false)
    lscene2 = LScene(lg[2,2], show_axis=false)
    plotmesh!(lscene2, mm;alpha=0.0, showsegments=true, segmentcolor=:darkgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene2, mm;color=Xs[2,:]./Ys[2,:], ceiling_offset=10, floor_offset=-15, colorrange=cr)
    viz!(lscene2, m_floor;color=Z)
    a,b = get_arrow(mm, directions[1])
    arrows3d!(lscene2, cmp-Point3f(a), Point3f(b))

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

function plot_directional_tuning!(lg, gidx::DirectionFiltered;kwargs...)
    # mild smoothing
    λ = get_direction_tuning(gidx;kwargs...)
    μr,ϕ = get_directional_tuning_strength(gidx;kwargs...)
    @show μr, ϕ
    pv = get_pvalue(gidx;kwargs...)
    for i in axes(λ,2)
        ax = PolarAxis(lg[i,1])
        hidedecorations!(ax)
        ax.thetagridvisible = true
        ax.rgridvisible = true
        ax.spinevisible = false
        ax.rticklabelsize = 10
        ax.thetaticklabelsize=10
        ax.rticksize = 0
        ax.thetaticksize = 0
        ax.thetaticklabelpad = 0
        lines!(ax, gidx.anglebins, λ[:,i], color=Cycled(i))
        ym = maximum(filter(isfinite, λ[:,i]))
        ym *= μr[i]
        linesegments!(ax,[ϕ[i],ϕ[i]], [0.0, ym], color=:red, linewidth=1)
        # indicate siginifance
        if pv[i] < 0.001
            lq = "**"
        elseif pv[i] < 0.01
            lq = "*"
        else
            lq = "ns"
        end
        ax.title = lq
        ax.titlesize = 12
        ax.titlegap = 0
        if i > 1
            rowgap!(lg, i-1, 0)
        end
        #Label(lg[i,2], lq, rotation=-π/2, tellheight=false,fontsize=14)
        #Label(lg[i,2], L"$μ_r = %$(round(μr[i], sigdigits=2))$ \\ $p = %$(round(pv[i], sigdigits=2))$", rotation=-π/2, tellheight=false, 
         #               fontsize=14)
    end
end

function plot_directional_tuning(gidx::DirectionFiltered;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(350,600))
        lg = GridLayout(fig[1,1])
        plot_directional_tuning!(lg, gidx;kwargs...)
        fig
    end
end