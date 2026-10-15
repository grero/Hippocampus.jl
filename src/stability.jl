using OptimalTransport

struct SpatialMapStability
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    ccm::Float64
    ccms::Vector{Float64}
end

function get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:all, smooth=false, α=0.5, niter=50,kwargs...)
    jm1 = JointMap(vpvrp, jocc,jocc_filtered;use_trials=use_trials)
    spm1 = Hippocampus.SpatialMapNew(jm1,m_floor)
    if smooth
        spml = SmoothedMap(spm1;method=:laplace, α=α, niter=niter)
        λ1 = get_rate_map(spml)
    else
        λ1 = get_rate_map(spm1)
    end
    λ1
end

function get_gaze_map(vpvrp, jocc, jocc_filtered, mm;use_trials=:all, smooth=false, α=0.5, niter=50,kwargs...)
    jm1 = JointMap(vpvrp, jocc,jocc_filtered;use_trials=use_trials)
    vm1 = Hippocampus.ViewMapNew(jm1,mm)
    if smooth
        vml = SmoothedMap(vm1;method=:laplace, α=α, niter=niter)
        λ1 = get_rate_map(vml)
    else
        λ1 = get_rate_map(vm1)
    end
    λ1
end

function get_map(vpvrp, jocc, jocc_filtered, mm;kwargs...)
    ed = embeddim(mm)
    if ed == 2
        return get_spatial_map(vpvrp, jocc, jocc_filtered,mm;kwargs...)
    elseif ed == 3
        return get_gaze_map(vpvrp, jocc, jocc_filtered,mm;kwargs...)
    end
    return nothing
end

function get_cross_correlation(λ1, λ2)
    qidx = (isfinite.(λ1)).&(isfinite.(λ2))
    cc = crosscor(λ1[qidx], λ2[qidx])
    cm = maximum(cc)
end

function compare_maps(λ1, λ2;normalize=true)
    fidx1 = isfinite.(λ1)
    fidx2 = isfinite.(λ2)
    fidx = fidx1.&fidx2
    if normalize
        λ1n =  (λ1.-mean(λ1[fidx1]))
        λ1n ./=norm(λ1n[fidx1])
        λ2n = (λ2.-mean(λ2[fidx2]))
        λ2n ./= norm(λ2n[fidx2])
    else
        λ1n = λ1
        λ2n = λ2
    end
    λ1n[fidx]'*λ2n[fidx]
end

function heat_kernel_cross_correlation(λ1, λ2, mm::SimpleMesh;σ=1.0)
    D = distancematrix(mm;between_centroids=false)
    heat_kernel_cross_correlation(λ1,λ2,D;σ=σ) 
end

function heat_kernel_cross_correlation(λ1, λ2, D::AbstractMatrix{<:Real};σ=1.0)
    total = 0.0
    μ1 = mean(filter(isfinite, λ1))
    μ2 = mean(filter(isfinite, λ2))
    for (k1,l1) in enumerate(λ1)
        if !isfinite(l1)
            continue
        end
        for (k2,l2) in enumerate(λ2)
            if !isfinite(l2)
                continue
            end
            d = D[k1,k2]
            kk = exp(-d^2/σ^2)
            total += (l1-μ1)*(l2-μ2)*kk
        end
    end
    total
end

function get_cross_correlation_1(λ1::AbstractVector{<:Real}, λ2::AbstractVector{<:Real}, mm::SimpleMesh;max_x_shift=nothing, max_y_shift=nothing)
    nn = nelements(mm)
    max_shift = div(round(Int64, ceil(sqrt(nn))),2)
    if max_x_shift === nothing
        max_x_shift = max_shift
    end
    if max_y_shift === nothing
        max_y_shift = max_shift
    end
    # do cross-correlation the brute force way, by explicitly shifting one map with respect to the other
    nn = nelements(mm)

    v = get_tangential_space(mm);
    A = adjacencymatrix(mm)
    λ1s = copy(λ1)
    qq = zeros(2*max_y_shift+1, 2*max_x_shift+1)
    shifts = Matrix{Tuple{Float64, Float64}}(undef, 2*max_y_shift+1, 2*max_x_shift+1)
    qq[max_y_shift+1, max_x_shift+1] = compare_maps(λ1, λ2) 
    shifts[max_y_shift+1, max_x_shift+1] = (0.0, 0.0)
    for (j,sx) in enumerate(1:1:max_x_shift)
        λ1s = shift_map(λ1, Vec2(-1.0, 0.0), v, A;nsteps=sx)
        qq[max_y_shift+1,j] = compare_maps(λ1s, λ2)
        shifts[max_y_shift+1,j] = (0.0, -sx)
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, -1.0), v, A;nsteps=sy)
            qq[i,j] = compare_maps(λ1ss, λ2)
            shifts[i,j] = (-sy,-sx)
        end
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, 1.0), v, A;nsteps=sy)
            qq[max_y_shift+i+1,j] = compare_maps(λ1ss, λ2)
            shifts[max_y_shift+i+1,j] = (-sx, sy)
        end
    end
     for (j,sx) in enumerate(1:1:max_x_shift)
        λ1s = shift_map(λ1, Vec2(1.0, 0.0), v, A;nsteps=sx)
        qq[max_y_shift+1,max_x_shift+j+1] = compare_maps(λ1s, λ2)
        shifts[max_y_shift+1,max_x_shift+j+1] = (0.0, sx)
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, -1.0), v, A;nsteps=sy)
            qq[i,max_x_shift+j+1] = compare_maps(λ1ss, λ2)
            shifts[i,max_x_shift+j+1] = (-sy,sx)
        end
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, 1.0), v, A;nsteps=sy)
            qq[max_y_shift+i+1,max_x_shift+j+1] = compare_maps(λ1ss, λ2)
            shifts[max_y_shift+i+1,max_x_shift+j+1] = (sy,sx)
        end
    end
    qq,shifts
end

function get_cross_correlation(λ1::AbstractVector{<:Real}, λ2::AbstractVector{<:Real}, mm::SimpleMesh;stepsize::Integer=1)
    # get the bins
    cp = Point2f.(Tuple.(centroid.(mm)))
    # get the binsize
    binsize = [v.val for v in sqrt.(measure.(mm))]

    D = maximum(norm.(cp .- permutedims(cp)))
    max_shift = round(Int64,ceil(floor(D/sqrt(2)/2)/binsize[1]))
    steps = -max_shift:stepsize:max_shift
    qq = zeros(length(steps), length(steps))
    lags = Matrix{Tuple{Float64, Float64}}(undef, size(qq)...)
    _binsize = binsize[1]
    for (i,sx) in enumerate(steps)
        λ1s,_ = shift_mass(λ1, cp, binsize, Vec2(sx*_binsize, 0.0)) 
        for (j,sy) in enumerate(steps)
            λ1ss,_ = shift_mass(λ1s, cp, binsize, Vec2(0.0, sy*_binsize))
            qq[j,i] = compare_maps(λ1ss, λ2) 
            lags[j,i] = (sy*_binsize,sx*_binsize)
        end
    end
    qq,lags
end

"""
    get_dissimilarity(λ1, λ2, mm::SimpleMesh)

Compute the dissimilarity between maps `λ1` and `λ2` defined on the mesh `mm`, using the sinkhole algorithm.
"""
function get_dissimilarity(λ1, λ2, mm::SimpleMesh)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    D = Hippocampus.distancematrix(mm;between_centroids=false)
    get_dissimilarity(λ1, λ2,D[fidx,fidx])
end

function get_dissimilarity(λ1, λ2, D::AbstractMatrix{<:Real})
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    f = λ1[fidx]./norm(λ1[fidx])
    g = λ2[fidx]./norm(λ2[fidx])
    P = sinkhorn_unbalanced(f,g,D, 1.0, 1.0, 0.01)
    #P = sinkhorn(f,g,D, 0.01)
    P,D
end

function get_null_distr(λ1, λ2, mm::SimpleMesh;nruns=1000)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    D = Hippocampus.distancematrix(mm;between_centroids=false,idx=fidx)
    λ1s = copy(λ1)
    λ2s = copy(λ2)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    dd = zeros(nruns)
    for r in 1:nruns
        for i in fidx
            if rand() < 0.5
                λ1s[i] = λ2[i]
            else
                λ1s[i] = λ1[i]
            end
            if rand() < 0.5
                λ2s[i] = λ1[i]
            else
                λ2s[i] = λ2[i]
            end
        end
        P,D = get_dissimilarity(λ1s, λ2s, D)
        dd[r] = sum(P.*D)
    end
    dd
end

function process_kwargs(::Type{SpatialMapStability},h::UInt32=zero(UInt32); nshuffles=1000, smooth=false, α=0.1, niter=50, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles => nshuffles),h)
    h = crc32c(string(:smooth=>smooth),h)
    if smooth
        h = crc32c(string(:α=>α),h)
        h = crc32c(string(:niter=>niter),h)
    end
    h
end

function SpatialMapStability(;redo=fname->false, do_save=true, nshuffles=1000, kwargs...)
    fname = "spatial_map_stability.jld2"
    h = process_kwargs(SpatialMapStability;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStability, fname)
    else
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        jocc,qdata = cd(DPHT.process_level("session")) do
            jocc = JointOccupancy(;kwargs...)
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
            jocc,qdata
        end
        jocc_filtered = JointFilteredOccupancy(jocc, qdata;kwargs...)

        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        obj = SpatialMapStability(vpvrp, jocc, jocc_filtered, m_floor;nshuffles=nshuffles, kwargs...)
        if do_save
            save_jld2(obj, fname;nshuffles=nshuffles, kwargs...)
        end
    end
    obj
end

struct SpatialMapStabilityCor
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cp1::Vector{Point2f}
    cp2::Vector{Point2f}
    cc::Float64
    ccs::Vector{Float64}
end

function process_kwargs(::Type{SpatialMapStabilityCor}, h::UInt32=zero(UInt32);nshuffles=1000, kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles))
    h
end

struct GazeMapStabilityCor
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cp1::Vector{Point2f}
    cp2::Vector{Point2f}
    cc::Float64
    ccs::Vector{Float64}
end

function process_kwargs(::Type{GazeMapStabilityCor}, h::UInt32=zero(UInt32);nshuffles=1000, kwargs...)
    h = process_kwargs(GazeResponseFields,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles))
    h
end


"""
    get_correspondence(rf1::SpatialResponseFields, rf2::SpatialResponseFields)

Compute maximum correspondence between the two spatial fields `rf1` and `rf2`
by finding the best alignment of each individual place field
"""
function get_correspondence(rf1::SpatialResponseFields, rf2::SpatialResponseFields)
    nrefinements = rf1.args[:nrefinements]
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    cp = Point2f.(Tuple.(centroid.(m_floor))) 
    binsize = [measure(m).val for m in m_floor]
    clusters1 = merge_fields(rf1)
    clusters2 = merge_fields(rf2)
    # attempt to align each field
    cp1 = Point2f[]
    cp2 = Point2f[]
    # What do we do if there are no clusters in either of these periods?
    for (i,cidx) in enumerate(clusters1)
        j = argmax(measure.(m_floor[rf1.binidx[cidx]]))
        push!(cp1, Point2(Tuple(centroid(m_floor[rf1.binidx[cidx[j]]]))))
    end
     for (i,cidx) in enumerate(clusters2)
        j = argmax(measure.(m_floor[rf2.binidx[cidx]]))
        push!(cp2, Point2(Tuple(centroid(m_floor[rf2.binidx[cidx[j]]]))))
    end

    ss = -Inf
    ij = (0,0)
    λ1 = rf1.λ
    λ2 = rf2.λ
    for (i,_cp1) in enumerate(cp1)
        for (j,_cp2) in enumerate(cp2)
            v = Vec2(_cp1 - _cp2)
            λ2s,_ = shift_mass(λ2, cp,binsize,  v)
            _ss = compare_maps(λ1, λ2s)
            if _ss > ss
                ss = _ss
                ij = (i,j)
            end
        end
    end
    ss, ij, cp1, cp2
end

struct SpatialMapStabilityGeo 
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cp1::Vector{Int64} # peak idx for map 1
    cp2::Vector{Int64} # peak idx for map 2
    cc::Float64
    ccs::Vector{Float64}
end

DPHT.filename(::Type{SpatialMapStabilityGeo}) = "spatial_map_stability_geo.jld2"

struct GazeMapStabilityGeo
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cp1::Vector{Int64}
    cp2::Vector{Int64}
    cc::Float64
    ccs::Vector{Float64}
end
DPHT.filename(::Type{GazeMapStabilityGeo}) = "gaze_map_stability_geo.jld2"

MapStabilityGeo = Union{SpatialMapStabilityGeo, GazeMapStabilityGeo}

get_value(X::T)  where T <: MapStabilityGeo = X.cc
get_surrogate_value(X::T) where T <: MapStabilityGeo = X.ccs
get_name(::Type{T}) where T <: MapStabilityGeo = "Geodesic distance"

function issignificant(X::T;pv_threshold=0.05) where T <: MapStabilityGeo
    ccs = get_surrogate_value(X)
    threshold = percentile(ccs, 100*pv_threshold)
    get_value(X) < threshold
end

function process_kwargs(::Type{GazeMapStabilityGeo}, h::UInt32=zero(UInt32);nshuffles=1000, kwargs...)
    h = process_kwargs(GazeResponseFields,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles))
    h
end

function process_kwargs(::Type{SpatialMapStabilityGeo}, h::UInt32=zero(UInt32);nshuffles=1000, kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles))
    h
end

get_response_field_type(::Type{GazeMapStabilityGeo}) = GazeResponseFields
get_response_field_type(::Type{SpatialMapStabilityGeo}) = SpatialResponseFields

struct GazeMapStabilityHK
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    c11::Float64 # self correlation
    c22::Float64 # self correlation
    c12::Float64 # cross-correlation
    c12ns::Vector{Float64} # normalized cross-correlation of surrogates
end

DPHT.filename(::Type{GazeMapStabilityHK}) = "gaze_maze_stability_hk.jld2"

struct SpatialMapStabilityHK
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    c11::Float64
    c22::Float64
    c12::Float64
    c12ns::Vector{Float64}
end
DPHT.filename(::Type{SpatialMapStabilityHK}) = "spatial_maze_stability_hk.jld2"

get_value(X::T)  where T <: MapStabilityHK = X.c12/sqrt(X.c11*X.c22)
get_surrogate_value(X::T) where T <: MapStabilityHK = X.c12ns
get_name(::Type{T}) where T <: MapStabilityHK = "Cross-correlation"

function issignificant(X::T;pv_threshold=0.05) where T <: MapStabilityHK
    ccs = get_surrogate_value(X)
    cc = get_value(X)
    threshold = percentile(ccs, 100*(1-pv_threshold))
    cc > threshold
end

get_response_field_type(::Type{GazeMapStabilityHK}) = GazeResponseFields
get_response_field_type(::Type{SpatialMapStabilityHK}) = SpatialResponseFields

MapStabilityHK = Union{GazeMapStabilityHK, SpatialMapStabilityHK}
MapStability = Union{MapStabilityHK, MapStabilityGeo}

GazeMapStabilityAll = Union{GazeMapStabilityGeo, GazeMapStabilityHK}
SpatialMapStabilityAll = Union{SpatialMapStabilityGeo, SpatialMapStabilityHK}

get_colormap(::Type{T}) where T <: SpatialMapStabilityAll = :rain
get_colormap(::Type{T}) where T <: GazeMapStabilityAll = :navia

function process_kwargs(::Type{T}, h::UInt32=zero(UInt32);nshuffles=1000, σ=1.0, kwargs...) where T<:MapStabilityHK
    Trf = get_response_field_type(T)
    h = process_kwargs(Trf,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    h = crc32c(string(:σ=>σ),h)
    h
end

function get_peak_idx(λ::AbstractVector{<:Real}, idx::Vector{T}) where T <: AbstractVector{T2} where T2 <: Integer
    c1 = fill(0, length(idx))
    for (ii,cidx1) in enumerate(idx)
        # find the elemnt associated with the larget value within cidx1
        c1[ii] = cidx1[argmax(λ[cidx1])]
    end
    c1
end

function get_average_geodesic_distance(rf1::T, rf2::T) where T <: AbstractResponseFields
     # identity the significant peaks in both maps
    nrefinements = rf1.args[:nrefinements]
    mm = get_mesh(T,nrefinements)
    clusters1 = Hippocampus.merge_fields(rf1)
    nclusters1 = Hippocampus.get_num_fields(rf1)
    cidx1 = findall(dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001)
    clusters1 = clusters1[cidx1]

    clusters2 = Hippocampus.merge_fields(rf2)
    nclusters2 = Hippocampus.get_num_fields(rf2)
    cidx2 = findall(dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001)
    clusters2 = clusters2[cidx2]

    c1 = get_peak_idx(rf1.λ, [rf1.binidx[c] for c in clusters1])
    c2 = get_peak_idx(rf2.λ, [rf2.binidx[c] for c in clusters2])
    get_average_geodesic_distance(mm, c1, c2)
end

function get_average_geodesic_distance(mm::SimpleMesh, peak1::AbstractVector{<:Integer}, peak2::AbstractVector{<:Integer})
    A = adjacencymatrix(mm)
    G = SimpleGraph(A)
    N = embeddim(mm)

    dd = zeros(length(peak2), length(peak1))
    for (j,_c1) in enumerate(peak1)
        dj = dijkstra_shortest_paths(G, _c1)
        for (i,_c2) in enumerate(peak2)
            pth = get_path(dj, _c2)
            dd[i,j] = Meshes.Unitful.ustrip(sum(norm.(diff(centroid.(mm[pth])))))
        end
    end
    pidx = argmin(dd, dims=1)
    dmin,pidx = dd[pidx], pidx
    dmin, pidx
end

function get_map_stability_geo(::Type{T};redo=fname->false, do_save=true, nshuffles=1000, kwargs...) where T <: MapStabilityGeo
    fname = DPHT.filename(T)
    h = process_kwargs(T;nshuffles=nshuffles, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
     if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilityCor, fname)
    else
        T_rf = get_response_field_type(T)
        rf1 = get_response_fields(T_rf, 1000;use_trials=:firstHalf, kwargs...)
        rf2 = get_response_fields(T_rf, 1000;use_trials=:secondHalf, kwargs...)
        nrefinements = rf1.args[:nrefinements]
        mm = get_mesh(T_rf,nrefinements)
        peaks1 = Hippocampus.get_peaks(rf1.λ, mm;t=2.5)
        c1 = Hippocampus.get_peak_idx(rf1.λ, peaks1)
        peaks2 = Hippocampus.get_peaks(rf2.λ, mm;t=2.5)
        c2 = Hippocampus.get_peak_idx(rf2.λ, peaks2)
        dmin,pidx = Hippocampus.get_average_geodesic_distance(mm, c1, c2)
        ss = mean(dmin)
        jocc,qdata,rpdata = cd(DPHT.process_level("session")) do
            jocc = JointOccupancy(;kwargs...)
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
            rp = RippleData()
            jocc,qdata,rp
        end
        jocc_filtered = JointFilteredOccupancy(jocc, qdata;kwargs...)

        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        rs2 = Hippocampus.RandomlyShiftedSpiketrains(;use_trials=:secondHalf, trial_start=2,redo=fname->false)

        smoothing_method = rf1.args[:smoothing_method]
        α = rf1.args[:α]
        niter = rf1.args[:niter]
        λ1 = rf1.λ
        ccs = zeros(nshuffles)
        trial_start = get(kwargs, :trial_start, 2)
        @showprogress for i in 1:nshuffles
            vpvrp = ViewAndPlaceRepresentationNew(rs2.timestamps[:,i]/1000.0, rpdata, qdata;trial_start=trial_start)
            λ2 = get_map(vpvrp, jocc, jocc_filtered, mm;use_trials=:secondHalf, smooth=true, α=α, niter=niter)
            peaks2 = Hippocampus.get_peaks(λ2, mm;t=2.5)
            _c2 = Hippocampus.get_peak_idx(λ2, peaks2)
            dmin,_ = Hippocampus.get_average_geodesic_distance(mm, c1, _c2)
            ccs[i] = mean(dmin) 
        end
        obj = T(λ1, rf2.λ, c1, c2, ss, ccs)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
        obj
    end
end

function get_map_stability_hk(::Type{T};redo=fname->false, do_save=true, nshuffles=1000, σ=1.0, kwargs...) where T <: MapStabilityHK
    fname = DPHT.filename(T)
    h = process_kwargs(T;nshuffles=nshuffles, σ=σ, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
     if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilityCor, fname)
    else
        T_rf = get_response_field_type(T)
        rf1 = get_response_fields(T_rf, 1000;use_trials=:firstHalf, kwargs...)
        rf2 = get_response_fields(T_rf, 1000;use_trials=:secondHalf, kwargs...)
        nrefinements = rf1.args[:nrefinements]
        mm = get_mesh(T_rf,nrefinements)
        D = distancematrix(mm;between_centroids=false)
        c11 = heat_kernel_cross_correlation(rf1.λ, rf1.λ, D;σ=σ)
        c22 = heat_kernel_cross_correlation(rf2.λ, rf2.λ, D;σ=σ)
        c12 = heat_kernel_cross_correlation(rf1.λ, rf2.λ, D;σ=σ)
        jocc,qdata,rpdata = cd(DPHT.process_level("session")) do
            jocc = JointOccupancy(;kwargs...)
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
            rp = RippleData()
            jocc,qdata,rp
        end
        jocc_filtered = JointFilteredOccupancy(jocc, qdata;kwargs...)

        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        rs2 = Hippocampus.RandomlyShiftedSpiketrains(;use_trials=:secondHalf, trial_start=2,redo=fname->false)

        smoothing_method = rf1.args[:smoothing_method]
        α = rf1.args[:α]
        niter = rf1.args[:niter]
        λ1 = rf1.λ
        ccs = zeros(nshuffles)
        trial_start = get(kwargs, :trial_start, 2)
        c12ns = zeros(nshuffles)
        @showprogress for i in 1:nshuffles
            vpvrp = ViewAndPlaceRepresentationNew(rs2.timestamps[:,i]/1000.0, rpdata, qdata;trial_start=trial_start)
            λ2 = get_map(vpvrp, jocc, jocc_filtered, mm;use_trials=:secondHalf, smooth=true, α=α, niter=niter)
            c22s = heat_kernel_cross_correlation(λ2, λ2,D;σ=σ)
            c12s = heat_kernel_cross_correlation(rf1.λ, λ2,D;σ=σ)
            c12ns[i] = c12s/sqrt(c11*c22s)
        end
        obj = T(rf1.λ, rf2.λ, c11, c22, c12, c12ns)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
        obj
    end
end

function get_map_stability(::Type{T}, args...;kwargs...) where T <: MapStabilityGeo
    get_map_stability_geo(T, args...;kwargs...)
end

function get_map_stability(::Type{T}, args...;kwargs...) where T <: MapStabilityHK
    get_map_stability_hk(T, args...;kwargs...)
end

function SpatialMapStabilityCor(;redo=fname->false, do_save=true, nshuffles=1000, kwargs...)
    fname = "spatial_map_stability_cor.jld2"
    h = process_kwargs(SpatialMapStabilityCor;nshuffles=nshuffles, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilityCor, fname)
    else
        rf1 = get_response_fields(SpatialResponseFields, 1000;use_trials=:firstHalf, kwargs...)
        rf2 = get_response_fields(SpatialResponseFields, 1000;use_trials=:secondHalf, kwargs...)
        ss,ij, cp1, cp2 = get_correspondence(rf1, rf2)
        v = Vec2(cp1[ij[1]] - cp2[ij[2]])
        # get significance by repeatedly scrambling one map, smoothing, then re-computing correspondenceo
        # using the original locations
        # why the original locations? Because we have no guarantee that the scrambled maps themselves have any
        # prominent local features
        # get the correspoding maps
        nrefinements = rf1.args[:nrefinements]
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        jocc,qdata,rpdata = cd(DPHT.process_level("session")) do
            jocc = JointOccupancy(;kwargs...)
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
            rp = RippleData()
            jocc,qdata,rp
        end
        jocc_filtered = JointFilteredOccupancy(jocc, qdata;kwargs...)

        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        rs2 = Hippocampus.RandomlyShiftedSpiketrains(;use_trials=:secondHalf, trial_start=2,redo=fname->false)

        smoothing_method = rf1.args[:smoothing_method]
        α = rf1.args[:α]
        niter = rf1.args[:niter]
         cp = Point2f.(Tuple.(centroid.(m_floor))) 
        binsize = [measure(m).val for m in m_floor]
        λ1 = rf1.λ
        ccs = zeros(nshuffles)
        trial_start = get(kwargs, :trial_start, 2)
        @showprogress for i in 1:nshuffles
            vpvrp = ViewAndPlaceRepresentationNew(rs2.timestamps[:,i]/1000.0, rpdata, qdata;trial_start=trial_start)
            λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:secondHalf, smooth=true, α=α, niter=niter)
            λ2s,_ = shift_mass(λ2, cp,binsize,  v)
            _ss = compare_maps(λ1, λ2s) 
            ccs[i] = _ss
        end
        obj = SpatialMapStabilityCor(λ1, rf2.λ, cp1, cp2, ss, ccs)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
        obj
    end
end

function get_spatial_map_stability(;nshuffles=10_000, kwargs...)
    # FIXME: There is something fishy with this function that prevents me from getting the result in the REPL
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    trial_start= get(kwargs, :trial_start, 2)
    # assume we are in the correct directory; load whatever we need
    sp = Spiketrain()

    rp = cd(DPHT.process_level(level(RippleData))) do
        RippleData()
    end
    rs1 = RandomlyShiftedSpiketrains(;nshifts=nshuffles, use_trials = :firstHalf, trial_start=trial_start)
    rs2 = RandomlyShiftedSpiketrains(;nshifts=nshuffles, use_trials = :secondHalf, trial_start=trial_start)

    unity_gaze_data = cd(DPHT.process_level("session")) do
        UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
    end
    jocc = cd(DPHT.process_level(level(JointOccupancy))) do
        JointOccupancy(;redo=fname->false, kwargs...)
    end
    jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
    vpvrp = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...) 

    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    P,D = get_dissimilarity(λ1, λ2, m_floor)
    ccm = sum(P.*D)
    ccms = zeros(nshuffles)
    for (ii,(sp1,sp2)) in enumerate(zip(eachcol(rs1.timestamps), eachcol(rs2.timestamps)))
        vpvrp1 = ViewAndPlaceRepresentationNew(sp1/1000.0,rp,unity_gaze_data;kwargs...) 
        λ1s = get_spatial_map(vpvrp1, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
        vpvrp2 = ViewAndPlaceRepresentationNew(sp2/1000.0,rp,unity_gaze_data;kwargs...) 
        λ2s = get_spatial_map(vpvrp2, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
        P,D = get_dissimilarity(λ1s, λ2s, m_floor)
        ccms[ii] = sum(P.*D)
    end
    ccm, ccms
end

struct SpatialMapStabilitySimple
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cc::Float64
    ccs::Vector{Float64}
end

function process_kwargs(::Type{SpatialMapStabilitySimple}, h::UInt32=zero(UInt32);nshuffles=10_000, stepsize::Integer=1,kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    if stepsize != 1
        h = crc32c(string(:stepsize=>stepsize),h)
    end
    h
end

function SpatialMapStabilitySimple(;redo=fname->false, do_save=true, nshuffles=10_000, stepsize=1, kwargs...)
    fname = "spatial_map_stability_simply.jld2"
    h = process_kwargs(SpatialMapStabilitySimple;nshuffles=nshuffles,kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilitySimple, fname)
    else
        rf1 = get_response_fields(SpatialResponseFields, 1000;use_trials=:firstHalf, pv_threshold=0.001, kwargs...)
        rf2 = get_response_fields(SpatialResponseFields, 1000;use_trials=:secondHalf, pv_threshold=0.001, kwargs...)
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        #jm1 =  JointMap(;use_trials=:firstHalf, kwargs...)
        #spm1 = SpatialMapNew(jm1,m_floor);
        #λ1 = get_rate_map(spm1);
        λ1 = rf1.λ
        #jm2 =  JointMap(;use_trials=:secondHalf, kwargs...)
        #spm2 = SpatialMapNew(jm2,m_floor); 
        #λ2 = get_rate_map(spm2);
        λ2 = rf2.λ

        qq, lags = get_cross_correlation(λ1, λ2, m_floor;stepsize=stepsize)
        ccq = maximum(qq)
        ccqs = zeros(nshuffles)
        @showprogress for ii in 1:nshuffles
            λ2s = shuffle(λ2)
            qq, _ = get_cross_correlation(λ1, λ2s, m_floor;stepsize=stepsize)
            ccqs[ii] = maximum(qq)
        end
        ccq, ccqs
        obj = SpatialMapStabilitySimple(λ1, λ2,ccq, ccqs)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
    end
    return obj
end


struct SpatialMapStabilityNew
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cc_ff::Vector{Float64}
    cc_fs::Vector{Float64}
    ds_ff::Vector{Float64}
    ds_fs::Vector{Float64}
end

function process_kwargs(::Type{SpatialMapStabilityNew}, h::UInt32=zero(UInt32);nshuffles=10_000, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    h
end


function SpatialMapStabilityNew(;redo=fname->false, do_save=true, nshuffles=10_000, kwargs...)
    fname = "spatial_map_stability_new.jld2"
    h = process_kwargs(SpatialMapStabilityNew;nshuffles=nshuffles, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilityNew, fname)
    else

        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))

        sp = Spiketrain()

        rp = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end

        unity_gaze_data = cd(DPHT.process_level("session")) do
            UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
        end
        jocc = cd(DPHT.process_level(level(JointOccupancy))) do
            JointOccupancy(;redo=fname->false, kwargs...)
        end
        nt = length(jocc.index)
        nthalf = div(nt,2)
        trialidx = [1:nt;]
        trialidx = collect(1:nt)
        jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
        vpvrp = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...) 

        cc_ff = zeros(nshuffles)
        ds_ff = zeros(nshuffles)
        cc_fs = zeros(nshuffles)
        ds_fs = zeros(nshuffles)
        @showprogress for i in 1:nshuffles
            # random maps from the first half of the trials
            trialidx1 = sort(shuffle(1:nthalf)[1:div(nthalf,2)])
            trialidx2 = setdiff(1:nthalf, trialidx1)
            λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx1,kwargs...)
            λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)
            P,D = get_dissimilarity(λ1, λ2, m_floor)
            ds_ff[i] = sum(P.*D)
            qq, lags = get_cross_correlation(λ1, λ2, m_floor)
            cc_ff[i] = maximum(qq)

            #second half
            trialidx2 = sort(shuffle(nthalf:1+nt)[1:div(nt-nthalf+1,2)])
            λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)

            P,D = get_dissimilarity(λ1, λ2, m_floor)
            ds_fs[i] = sum(P.*D)
            qq, lags = get_cross_correlation(λ1, λ2, m_floor)
            cc_fs[i] = maximum(qq)
        end
        λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
        λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
        obj = SpatialMapStabilityNew(λ1, λ2,cc_ff, cc_fs, ds_ff, ds_fs)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
    end
    obj
end

function get_spatial_map_stability(vpvrp, jocc, jocc_filtered, m_floor::SimpleMesh;kwargs...)
    # compute stability of mean spatial response for 1st and 2nd half of the trials, compare with equivalent number of random trials
    # If there is a systematic shift 
    nt = length(jocc.index)
    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    P,D = get_dissimilarity(λ1, λ2, m_floor)
    ccm = sum(P.D)
end

function SpatialMapStability(vpvrp, jocc, jocc_filtered, m_floor::SimpleMesh;nshuffles=1000, kwargs...)
    # compute stability of mean spatial response for 1st and 2nd half of the trials, compare with equivalent number of random trials
    # If there is a systematic shift 
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    nt = length(jocc.index)
    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    qq,lags = get_cross_correlation(λ1, λ2, m_floor)
    ccm = maximum(qq)
    ccms = zeros(nshuffles)
    nthalf = div(nt,2)
    trialidx = [1:nt;]
    @showprogress for i in 1:length(ccms)
        #trialidx1 = sort(shuffle(trialidx)[1:nthalf])
        #trialidx2 = setdiff(trialidx, trialidx1)
        #λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx1,kwargs...)
        #λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)
        qqs,lags = get_cross_correlation(shuffle(λ1), shuffle(λ2),m_floor)
        ccms[i] = maximum(qqs)
    end
    SpatialMapStability(λ1, λ2, ccm, ccms)
end

## plots

function plot_stability_summary(::Type{SpatialMapStabilitySimple}, celldirs::Vector{String};_plot_theme=plot_theme, kwargs...)
    kargs = (nshuffles=1000, nrefinements=(p=2,g=2), min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, smooth=false)
    # TODO: Also add the spatial maps here
    spm_stability_simple = map(celldirs) do celldir
        spm_stability = cd(celldir) do
            spm = Hippocampus.SpatialMapStabilitySimple(;kargs...)
            spm
        end
        spm_stability.cc, spm_stability.ccs
    end
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=2))
    spm_stability = [s[1] for s in spm_stability_simple]
    spm_stability_s = [percentile(s[2], 95) for s in spm_stability_simple]
    vidx = spm_stability .> spm_stability_s
    midx = findall(vidx)
    @show length(midx)
    idx0 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 5)))]
    idx1 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 50)))]
    idx2 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 95)))]
    @show idx0 idx1 idx2

    with_theme(_plot_theme) do
        fig = Figure(size=(1.5*650,1.5*300))
        ax = Axis(fig[1,1])
        hist!(ax, spm_stability[vidx.==false];color=(:white, 0.0), strokecolor=:gray, strokewidth=2.0)
        hist!(ax, spm_stability[vidx.==true];color=(:white,0.0), strokecolor=:royalblue4, strokewidth=2.0)
        ax.xlabel = "Stability"
        # show example of the 5th percentile, the median, and the 95th percentile cells in terms of stability
        lg = GridLayout(fig[1,2])
        lg1 = GridLayout(lg[1,1], alignmode=Outside(5))
        lg2 = GridLayout(lg[1,2], alignmode=Outside(5))
        lg3 = GridLayout(lg[1,3], alignmode=Outside(5))
        # indicate these points on the histogram
        colors = [:pink, :red, :orange]
        vlines!(ax, spm_stability[[idx0,idx1,idx2]], color=colors)
        # draw boxes around the corresponding plots
        Makie.Box(lg[1,1], color=(:white, 0.0), strokecolor=colors[1])
        Makie.Box(lg[1,2], color=(:white,0.0), strokecolor=colors[2])
        Makie.Box(lg[1,3], color=(:white,0.0), strokecolor=colors[3])
        for (ii,(idx,_lg)) in enumerate(zip([idx0, idx1, idx2],[lg1,lg2,lg3]))
            jm1,jm2 = cd(celldirs[idx])  do
                jm1 = JointMap(;use_trials=:firstHalf,kargs...)
                jm2 = JointMap(;use_trials=:secondHalf, kargs...)
                jm1, jm2
            end
            spm1 = SpatialMapNew(jm1, m_floor)
            spm2 = SpatialMapNew(jm2, m_floor)
            ax1 = Axis(_lg[1,1], aspect=1.0)
            ax2 = Axis(_lg[2,1], aspect=1.0)
            for _ax in [ax1, ax2]
                hidedecorations!(_ax)
                _ax.bottomspinevisible = false
                _ax.leftspinevisible = false
            end
            λ1 = get_rate_map(spm1)
            λ2 = get_rate_map(spm2)
            cr = extrema(filter(isfinite, [λ1;λ2]))
            viz!(ax1, m_floor;color=:darkgray)
            viz!(ax1, m_floor;color=get_rate_map(spm1),colormap=:rain, colorrange=cr)
            viz!(ax2, m_floor;color=:darkgray)
            viz!(ax2, m_floor;color=get_rate_map(spm2),colormap=:rain, colorrange=cr)

            rf = cd(celldirs[idx]) do
                rf = get_response_fields(Hippocampus.SpatialResponseFields, 1000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001, use_trials=:all)
            end
            bb = find_boundaries(rf)
            # superimpose the boundaries from the original place cell calculation for illustration purposes only
            for _bb in bb
                viz!(ax1, _bb, color=:red)
                viz!(ax2, _bb, color=:red)
            end
           
        end
        lgf = GridLayout(fig[1,3])
        Label(lgf[1,1], "First half", rotation=π/2, tellheight=false)
        Label(lgf[2,1], "Second half", rotation=π/2, tellheight=false)
        colsize!(fig.layout, 1, Relative(0.4))
        fig
    end
end

function plot_stability_summary(::Type{SpatialMapStabilityCor}, celldirs::Vector{String};_plot_theme=plot_theme, kwargs...)
    kargs = (nshuffles=1000, nrefinements=(p=3,g=2), min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50)
    # TODO: Also add the spatial maps here
    spm_stability_simple = map(celldirs) do celldir
        spm_stability = cd(celldir) do
            spm = Hippocampus.SpatialMapStabilityCor(;kargs...)
            spm
        end
        spm_stability.cc, spm_stability.ccs
    end
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=3))
    spm_stability = [s[1] for s in spm_stability_simple]
    spm_stability_s = [percentile(s[2], 95) for s in spm_stability_simple]
    vidx = spm_stability .> spm_stability_s
    midx = findall(vidx)
    @show length(midx)
    idx0 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 5)))]
    idx1 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 50)))]
    idx2 = midx[argmin(norm.(spm_stability[vidx] .- percentile(spm_stability[vidx], 95)))]
    @show idx0 idx1 idx2

    with_theme(_plot_theme) do
        fig = Figure(size=(1.5*650,1.5*300))
        ax = Axis(fig[1,1])
        hist!(ax, spm_stability[vidx.==false];color=(:white, 0.0), strokecolor=:gray, strokewidth=2.0)
        hist!(ax, spm_stability[vidx.==true];color=(:white,0.0), strokecolor=:royalblue4, strokewidth=2.0)
        ax.xlabel = "Stability"
        # show example of the 5th percentile, the median, and the 95th percentile cells in terms of stability
        lg = GridLayout(fig[1,2])
        lg1 = GridLayout(lg[1,1], alignmode=Outside(5))
        lg2 = GridLayout(lg[1,2], alignmode=Outside(5))
        lg3 = GridLayout(lg[1,3], alignmode=Outside(5))
        # indicate these points on the histogram
        colors = [:pink, :red, :orange]
        vlines!(ax, spm_stability[[idx0,idx1,idx2]], color=colors)
        # draw boxes around the corresponding plots
        Makie.Box(lg[1,1], color=(:white, 0.0), strokecolor=colors[1])
        Makie.Box(lg[1,2], color=(:white,0.0), strokecolor=colors[2])
        Makie.Box(lg[1,3], color=(:white,0.0), strokecolor=colors[3])
        for (ii,(idx,_lg)) in enumerate(zip([idx0, idx1, idx2],[lg1,lg2,lg3]))
            rf1,rf2 = cd(celldirs[idx])  do
                rf1 = get_response_fields(SpatialResponseFields, 1000;use_trials=:firstHalf,kargs...)
                rf2 = get_response_fields(SpatialResponseFields, 1000;use_trials=:secondHalf,kargs...)
                rf1, rf2
            end
            _lg1 = GridLayout(_lg[1,1])
            _lg2 = GridLayout(_lg[2,1])
            cr = extrema(filter(isfinite, [rf1.λ;rf2.λ]))
            ax1 = plot_response_fields!(_lg1, rf1;colorrange=cr, colormap=:rain, show_colorbar=false) 
            ax2 =plot_response_fields!(_lg2, rf2;colorrange=cr, colormap=:rain, show_colorbar=false) 

            rf = cd(celldirs[idx]) do
                rf = get_response_fields(Hippocampus.SpatialResponseFields, 1000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001, use_trials=:all)
            end
            bb = find_boundaries(rf)
            # superimpose the boundaries from the original place cell calculation for illustration purposes only
            for _bb in bb
                viz!(ax1, _bb, color=:red)
                viz!(ax2, _bb, color=:red)
            end
           
        end
        lgf = GridLayout(fig[1,3])
        Label(lgf[1,1], "First half", rotation=π/2, tellheight=false)
        Label(lgf[2,1], "Second half", rotation=π/2, tellheight=false)
        colsize!(fig.layout, 1, Relative(0.4))
        fig
    end
end

function plot_stability(args...;_plot_theme=plot_theme,kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(600,400))
        lg = GridLayout(fig[1,1])
        plot_stability!(lg, args...;kwargs...)
        fig
    end
end

function plot_stability!(lg, vpvrp::ViewAndPlaceRepresentationNew, jocc, jocc_filtered, m_floor::SimpleMesh;kwargs...)
    λ1 = Hippocampus.get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:firstHalf)
    λ2 = Hippocampus.get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:secondHalf)
    plot_stability!(lg, λ1, λ2, m_floor;kwargs...)
end

function plot_stability!(lg, spm::SpatialMapStabilityNew, m_floor::SimpleMesh;plot_cross_correlation=true)
    if plot_cross_correlation
        (s1,s2) = (spm.cc_ff, spm.cc_fs)
    else
        (s1,s2) = (spm.ds_ff, spm.ds_fs)
    end
    plot_stability!(lg, spm.λ1, spm.λ2, m_floor, s1,s2)
end

function plot_stability!(lg, spm::Union{SpatialMapStabilitySimple,SpatialMapStabilityCor}, m_floor::SimpleMesh)
    plot_stability!(lg, spm.λ1, spm.λ2, m_floor)
    Label(lg[1,1,Top()], "First half")
    Label(lg[1,2,Top()], "Second half")
    ax = Axis(lg[2,1:3])
    boxplot!(ax, fill(1.0, length(spm.ccs)), spm.ccs, color=:gray, orientation=:horizontal)
    vlines!(ax, spm.cc, linestyle=:dot)
    ax.xlabel = "1st vs 2nd half corr"
    ax.leftspinevisible = false
    ax.yticksvisible = false
    ax.yticklabelsvisible = false
    rowsize!(lg, 2, Relative(0.1))
end

function plot_stability!(lg, λ1::AbstractVector{<:Real},λ2::AbstractVector{<:Real}, m_floor, s1::AbstractVector{<:Real}, s2::AbstractVector{<:Real};ylabel="Cross-correlation", kwargs...)
    lg1 = GridLayout(lg[1,1])
    plot_stability!(lg1, λ1, λ2, m_floor;kwargs...)
    ax = Axis(lg1[2,1:2])
    Label(lg1[1,1,Top()], "First half")
    Label(lg1[1,2,Top()], "Second half")
    xx = [fill(1.0, length(s1));fill(2.0, length(s2))]
    yy = [s1;s2]
    fidx = isfinite.(yy)
    rainclouds!(ax, xx[fidx], yy[fidx])
    ax.ylabel = ylabel
    ax.xticklabelsvisible = true 
    ax.xticksvisible =  false 
    ax.bottomspinevisible = false 
    ax.xticks = ([1,2], ["fh→fh","fh→sh"])
    Label(lg1[1,1,TopLeft()], "A")
    Label(lg1[2,1,TopLeft()], "B")
end

function plot_stability!(lg, λ1::AbstractVector{<:Real},λ2::AbstractVector{<:Real}, m_floor;kwargs...)
    cr = extrema(filter(isfinite, [λ1;λ2]))
    colormap = get(kwargs, :colormap, :rain)
    ax1 = Axis(lg[1,1], aspect=1.0)
    viz!(ax1, m_floor;color=:darkgray)
    viz!(ax1, m_floor;color=λ1, colormap=colormap, colorrange=cr)

    ax2 = Axis(lg[1,2], aspect=1.0)
    viz!(ax2, m_floor;color=:darkgray)
    viz!(ax2, m_floor;color=λ2, colormap=colormap, colorrange=cr)

    Colorbar(lg[1,3], colorrange=cr, colormap=colormap, label="Firing rate [Hz]")
    for _ax in [ax1, ax2]
        hidedecorations!(_ax)
        _ax.bottomspinevisible = false
        _ax.leftspinevisible = false
    end
end

function plot_stability(celldir::String;kargs...)
    
end

function plot_stability!(lg, stab::GazeMapStabilityGeo, mm;kwargs...)
    # because of a bug in how the cp2 was saved, we nned to recompute this
    peaks2 = get_peaks(stab.λ2, mm;t=2.5)
    cp2 = get_peak_idx(stab.λ2, peaks2)
    # match the peaks by plotting them in the same color
    _, pidx = get_average_geodesic_distance(mm, stab.cp1, cp2)
    colors = [:firebrick1, :orange, :sienna, :hotpink3]
    @debug pidx
    @debug length(stab.cp1) length(cp2)
    @debug getindex.(pidx,1)
    @debug getindex.(pidx,2)
    if length(cp2) < length(stab.cp1)
        colors1 = vec(colors[getindex.(pidx,2)])
        colors2 = colors[1:length(cp2)]
    else
        colors1 = vec(colors[1:length(stab.cp1)])
        colors2 = fill(:gray45, length(cp2))
        colors2[getindex.(pidx,1)] .= colors1
    end
    # TODO: indicate geodesics
    colormap = get(kwargs, :colormap, :rain)
    lscene1 = LScene(lg[1,1],show_axis=false)
    Label(lg[1,1,Top()], "First half", tellwidth=false)
    plotmesh!(lscene1, mm;color=stab.λ1, floor_offset=0, ceiling_offset=15, showsegments=true, colormap=colormap)
    cr = extrema(filter(isfinite,stab.λ1))
    Colorbar(lg[1,2], colormap=colormap, colorrange=cr, label = "Firing rate [Hz]")
    viz!(lscene1, centroid.(mm[stab.cp1]), pointsize=10, color=colors1)
    lscene2 = LScene(lg[1,3],show_axis=false)
    Label(lg[1,3,Top()], "Second half", tellwidth=false)
    plotmesh!(lscene2, mm;color=stab.λ2, floor_offset=0, ceiling_offset=15, showsegments=true, colormap=colormap)
    cr = extrema(filter(isfinite,stab.λ2))
    @show colors2
    viz!(lscene2, centroid.(mm[cp2]),pointsize=10, color=colors2)
    Colorbar(lg[1,4], colormap=colormap, colorrange=cr, label="Firing rate [Hz]")
    # indicate distrubion
    ax3 = Axis(lg[1,5])
    boxplot!(ax3, fill(1.0, length(stab.ccs)), stab.ccs;color=:gray45)
    hlines!(ax3, stab.cc, color=:royalblue4, linestyle=:dot, linewidth=2.0)
    colsize!(lg, 5, 50)
    ax3.ylabel = "Geodesic dist [unit]"
    ax3.xticklabelsvisible = false
    ax3.xticksvisible = false
    ax3.bottomspinevisible = false
end

function plot_stability!(lg, stab::GazeMapStabilityGeo, mm, rf::GazeResponseFields;kwargs...)
    lg0 = GridLayout(lg[1,0])
    Label(lg[1,0,Top()], "Whole session", tellwidth=false)
    plot_response_fields!(lg0,rf;floor_offset=0.0, mazecolor=nothing, kwargs...)
    plot_stability!(lg, stab,mm;kwargs...)
end

function plot_stability(stab::GazeMapStabilityGeo, mm, args...;_plot_theme=plot_theme, kwargs...)
    with_theme(_plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_stability!(lg, stab, mm, args...;kwargs...)
        fig
    end
end

function plot_stability_summary!(lg, ::Type{T}, celldirs::Vector{String};kwargs...) where T <: MapStability
    kargs = (nshuffles=1000, nrefinements=(p=3,g=2), min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, smooth=true,smoothing_method=:laplace, α=0.1, niter=50)
    nv = length(celldirs)
    vm_stability_value = zeros(nv)
    vm_stability_sig = fill(false, nv)

    for (ii, celldir) in enumerate(celldirs)
       vm_stability = cd(celldir) do
       Hippocampus.get_map_stability(T;kargs...)
       end
       vm_stability_value[ii] = get_value(vm_stability)
       vm_stability_sig[ii] = issignificant(vm_stability)
    end

    ax = Axis(lg[1,1,])
    vms_sig = vm_stability_value[vm_stability_sig]
    @show sum(vm_stability_sig)
    l,m,h = percentile(vms_sig, [10,50,90])
    x0,idx0 = findmin(x->norm(x-l), vms_sig)
    x1,idx1 = findmin(x->norm(x-m), vms_sig)
    x2,idx2 = findmin(x->norm(x-h), vms_sig)
    idx0,idx1,idx2 = findall(vm_stability_sig)[[idx0,idx1,idx2]]
    @show idx0 idx1 idx2
    h1 = hist!(ax, vm_stability_value[(!).(vm_stability_sig)], color=:gray25)
    h2 = hist!(ax, vm_stability_value[vm_stability_sig], color=:royalblue4)
    # fit density, most just for visualls
    # scale density by height
    max1 = maximum([v[2] for v in h1.points.value[]])
    max2 = maximum([v[2] for v in h2.points.value[]])
    kde1 = kde(vm_stability_value[(!).(vm_stability_sig)])
    kde2 = kde(vm_stability_value[vm_stability_sig])
    lines!(ax, kde1.x, max1*kde1.density/maximum(kde1.density), color=:gray45)
    lines!(ax, kde2.x, max2*kde2.density/maximum(kde2.density), color=:royalblue2)
    colors = [:pink, :red, :orange]
    vlines!(ax, vm_stability_value[[idx0,idx1,idx2]], color=colors)
    ax.xlabel = get_name(T)
    ax.ylabel = "Count"

    lgq = GridLayout(lg[1,2])
    lg1 = GridLayout(lgq[1,1], alignmode=Outside(5))
    lg2 = GridLayout(lgq[1,2], alignmode=Outside(5))
    lg3 = GridLayout(lgq[1,3], alignmode=Outside(5))
    # indicate these points on the histogram
    # draw boxes around the corresponding plots
    Makie.Box(lgq[1,1], color=(:white, 0.0), strokecolor=colors[1])
    Makie.Box(lgq[1,2], color=(:white,0.0), strokecolor=colors[2])
    Makie.Box(lgq[1,3], color=(:white,0.0), strokecolor=colors[3])
    Trf = get_response_field_type(T)
    colormap = get_colormap(T)
    for (ii,(idx,_lg)) in enumerate(zip([idx0, idx1, idx2],[lg1,lg2,lg3]))
        rf1,rf2 = cd(celldirs[idx])  do
            rf1 = get_response_fields(Trf, 1000;use_trials=:firstHalf,kargs...)
            rf2 = get_response_fields(Trf, 1000;use_trials=:secondHalf,kargs...)
            rf1, rf2
        end
        _lg1 = GridLayout(_lg[1,1])
        _lg2 = GridLayout(_lg[2,1])
        cr = extrema(filter(isfinite, [rf1.λ;rf2.λ]))
        ax1 = plot_response_fields!(_lg1, rf1;colorrange=cr, colormap=colormap, show_colorbar=false, mazecolor=nothing,floor_offset=-20) 
        ax2 =plot_response_fields!(_lg2, rf2;colorrange=cr, colormap=colormap, show_colorbar=false, mazecolor=nothing, floor_offset=-20) 

        rf = cd(celldirs[idx]) do
            rf = get_response_fields(Trf, 1000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001, use_trials=:all)
        end
        bb = find_boundaries(rf)
        # superimpose the boundaries from the original place cell calculation for illustration purposes only
        for _bb in bb
            viz!(ax1, _bb, color=:red)
            viz!(ax2, _bb, color=:red)
        end
        
    end
    lgf = GridLayout(lg[1,3])
    Label(lgf[1,1], "First half", rotation=π/2, tellheight=false)
    Label(lgf[2,1], "Second half", rotation=π/2, tellheight=false)
    colsize!(lg, 1, Relative(0.3))
end

function plot_stability_summary(::Type{T}, celldirs::Vector{String};_plot_theme=plot_theme, kwargs...) where T <: MapStability
    w = 1024
    h = (7/16)*w
    with_theme(_plot_theme) do
        fig = Figure(size=(w,h))
        lg = GridLayout(fig[1,1])
        plot_stability_summary!(lg,T, celldirs;kwargs...)
        fig
    end
end