find_place_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_place_selective, celldirs;kwargs...)
find_view_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_view_selective, celldirs;kwargs...)

function find_selective_cells(func::Function, celldirs::Vector{String};skip_error=true, kwargs...)
    _is_selective = fill(false, length(celldirs))
    for (i,celldir) in enumerate(celldirs)
        try
            _is_selective[i] = func(celldir;kwargs...)
        catch ee
            if !skip_error
                rethrow(ee)
            end
        finally
            continue
        end
    end
    _is_selective
end

function is_place_selective(celldir::String;filtering="All", raycast_type="1px")
    fname=joinpath(celldir, "Filt$(filtering)",raycast_type,"vmpc.mat")
    if !ispath(fname)
        error("No vmpc object found at $fname")
    end

    vp = MAT.matread(fname)
    vp_data = vp["vmp"]["data"]
    is_selective(vp_data)
end

function is_view_selective(celldir::String;filtering="All", raycast_type="1px")
    fname=joinpath(celldir, "Filt$(filtering)",raycast_type,"vmsv.mat")
    if !ispath(fname)
        error("No vmsv object found at $fname")
    end

    vp = MAT.matread(fname)
    vp_data = vp["vms"]["data"]
    is_selective(vp_data)
end
    
function is_selective(vp_data::Dict{String,Any})
    # check for 
   sic_value, sic_threshold = get_sic_with_shuffle(vp_data) 
    return sic_value > sic_threshold 
end

function get_sic_with_shuffle(vp_data::Dict{String, Any})
    if "SIC" in keys(vp_data)
        sic_value = vp_data["SIC"]
    elseif "SIC_adsm" in keys(vp_data)
        sic_value = vp_data["SIC_adsm"]
    else
        error("No SIC variable found")
    end
    if "SICsh" in keys(vp_data)
        sic_threshold = percentile(vec(vp_data["SICsh"]), 95)
    elseif "crit_sm" in keys(vp_data)
        sic_threshold = vp_data["crit_sm"]
    else
        error("No shuffled data found")
    end
    sic_value, sic_threshold
end

function get_sic_with_shuffle(celldir::String;filtering="All", raycast_type="1px")
    fname=joinpath(celldir, "Filt$(filtering)",raycast_type,"vmpc.mat")
    if !ispath(fname)
        error("No vmpc object found at $fname")
    end

    vp = MAT.matread(fname)
    vp_data = vp["vmp"]["data"]
    sic_value, sic_threshold = get_sic_with_shuffle(vp_data)
    sic_value, sic_threshold
end

function plot_place_fields_with_outline(sessionnr=16)
    place_selective_cells = open("data/place_selective_cells.txt") do fid
       readlines(fid)
    end
    sessiondirs = Hippocampus.DPHT.get_level_path.("session", place_selective_cells)
    sessiondir = sessiondirs[sessionnr]
    cidx = findall(sessiondirs.==sessiondir)
    f = get_spatial_maps(place_selective_cells[cidx])

    # get the points of the place fields
    patches = [Hippocampus.field_outline(f[:,:,i];t=1.65) for i in 1:size(f,3)]
    # create outlines via convex hull
    xbins = range(-12.5f0, stop=12.5f0, length=40);
    ybins = range(-12.5f0, stop=12.5f0, length=40);
    all_hulls = Any[]
    for patch in patches
        _hulls = Any[]
        for p in patch
            points = [Meshes.Point(xbins[ci.I[1]], ybins[ci.I[2]]) for ci in  p]
            chull = convexhull(points)
            push!(_hulls, chull)
        end
        push!(all_hulls, _hulls)
    end
    # find points for with overlaps
    udata = cd(sessiondir) do
        UnityData()
    end
    spa = map(place_selective_cells[cidx]) do celldir
       cd(celldir) do
        Hippocampus.TrialAlignedSpiketrain()
       end
    end
    weight,bins = Hippocampus.compute_psth(spa, 0.05,8;tmax=20.0)
    hulls = cat(all_hulls..., dims=1)
    timestamps,pos = get_timepoints(udata, hulls)
    idx = CartesianIndex[]
    # find the bin index correpsonding to these time stamps
    allpos = Tuple{Float64,Float64}[]
    for i in 1:nt
        for (pp,tt) in zip(pos[i], timestamps[i])
            kk = searchsortedfirst(bins[1:size(weight,1)], tt)
            if 0 < kk <= size(weight,1)
                push!(idx, CartesianIndex(kk,i))
                push!(allpos, pp)
            end
        end
    end

    XX = zeros(length(idx), size(weight,3))
    for (j,ii) in enumerate(idx)
       XX[j,:] = weight[ii.I[1], ii.I[2],:]
    end

    fig = Figure()
    ax = Axis(fig[1,1])
    m_floor = Hippocampus.floor_topology3()
    colors = to_colormap(:tab10)
    viz!(ax, m_floor,color=:lightgray)
    for (ci,hull) in enumerate(all_hulls)
        for _hull in hull
            viz!(ax,_hull,color=colors[ci])
        end
    end
    fig
end

function plot_trajectories(udata::UnityData)
    nt = numtrials(udata)
    trajectories = Dict{Tuple{Int64, Int64}, Vector{Vector{Tuple{Float64,Float64}}}}()
    prev_posterid = 0
    for i in 1:nt
        # make sure the trial was correct
        if !(30 < udata.triggers[i,3] < 40)
            posterid = 0
        else
            posterid = udata.triggers[i,1] - 10
            kk = (prev_posterid, posterid)
            if !(kk in keys(trajectories))
                trajectories[kk] = Vector{Tuple{Float64, Float64}}[]
            end
            tg,posx,posy,hd = get_trial(udata, i;trial_start=2)
            push!(trajectories[kk], [(px,py) for (px,py) in zip(posx,posy)])
        end
        prev_posterid = posterid
    end
    m_floor = floor_topology3() 
    with_theme(plot_theme) do 
        fig = Figure()
        axes = [Axis(fig[i,j]) for i in 1:6, j in 1:6]
        hidedecorations!.(axes)
        _keys = collect(keys(trajectories))
        sort!(_keys)
        for k in _keys
            if (0 in k) || (k[1]==k[2])
                continue
            end
            ax = axes[k[1], k[2]]
            viz!(ax, m_floor, color=:lightgray)
            for v in trajectories[k]
                lines!(ax, Point2f.(v),color=:black)
            end
            scatter!(ax, Point2f(first(trajectories[k])[1]), color=:green)
        end
        fig
    end
end