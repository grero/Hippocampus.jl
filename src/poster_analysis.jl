poster_theme = Theme(Axis=(xlabelsize=36, ylabelsize=36,
                           xticklabelsize=36, yticklabelsize=36,
                           topspinevisible=false, rightspinevisible=false,
                           xgridvisible=false, ygridvisible=false,ylabelvisible=true,
                           xticklabelsvisible=true, xlabelvisible=true),
                     Scatter=(markersize=10px,),
                     Lines=(linewidth=3,),
                     fontsize=36)



find_place_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_place_selective, celldirs;kwargs...)
find_view_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_view_selective, celldirs;kwargs...)

function get_maze_category_colors()
    cmap = to_colormap(:tab20)
    _colors = reshape(cmap[[1:4;9:10;17:18;]], 2, 4)
    pillar_colors = HSV{Float32}[]
    for i in axes(_colors,2)
        append!(pillar_colors, range(HSV(_colors[1,i]), stop=HSV(_colors[2,i]), length=4))
    end
    floor_color = cmap[5]
    ceiling_color = cmap[19]
    wall_colors = range(HSV(cmap[11]),stop=HSV(cmap[12]), length=4)
    [ceiling_color;floor_color;wall_colors;pillar_colors]
end

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

function plot_poster_decoding_results(;cell_examples=(poster_selective=176, previous_poster_selective=292, both=9))
    allcelldirs = open("/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/cell_list.txt") do fid
        readlines(fid)
    end
    data = JLD2.load("data/poster_id_population_decoding_results.jld2")
    perf_c = data["performance_current"]
    μ_perf_c = dropdims(mean(perf_c,dims=1),dims=1)
    perf_p = data["performance_previous"]
    μ_perf_p = dropdims(mean(perf_p,dims=1),dims=1)
    ncells = data["ncells_per_session"]
    nt = data["correct_trials_per_session"]
    contrib_current = data["cell_contrib_current"]
    contrib_previous = data["cell_contrib_previous"]
    current_poster_selectivity = data["current_poster_selectivity_strength"]
    previous_poster_selectivity = data["previous_poster_selectivity_strength"]

    with_theme(poster_theme) do
        width = 11.0*2.5*72
        height = width
        fig = Figure(size=(width, height))
        # indivvidual cell responses
        lg2 = GridLayout(fig[1,1])

        lg21 = GridLayout(lg2[1,1])
        Label(lg21[1,1,TopLeft()], "A")
        plot_raster_and_psth!(lg21, allcelldirs[cell_examples.previous_poster_selective];previous=true,xlabelvisible=false, xticklabelsvisible=false,binsize=0.1)
        rowsize!(lg21, 1, Relative(0.6))
        lg31 = GridLayout(lg2[2,1])
        rowsize!(lg31, 1, Relative(0.6))
        plot_raster_and_psth!(lg31, allcelldirs[cell_examples.previous_poster_selective];previous=false, binsize=0.1)

        lg22 = GridLayout(lg2[1,2])
        Label(lg22[1,1,TopLeft()], "B")
        rowsize!(lg22, 1, Relative(0.6))
        lg32 = GridLayout(lg2[2,2])
        rowsize!(lg32, 1, Relative(0.6))
        lg23 = GridLayout(lg2[1,3])
        Label(lg23[1,1,TopLeft()], "C")
        rowsize!(lg23, 1, Relative(0.6))
        lg33 = GridLayout(lg2[2,3])
        rowsize!(lg33, 1, Relative(0.6))
        plot_raster_and_psth!(lg22, allcelldirs[cell_examples.poster_selective];previous=true,ylabelvisible=false,xlabelvisible=false, xticklabelsvisible=false,binsize=0.1)
        plot_raster_and_psth!(lg32, allcelldirs[cell_examples.poster_selective];previous=false, ylabelvisible=false, binsize=0.1)

        plot_raster_and_psth!(lg23, allcelldirs[cell_examples.both];previous=true, ylabelvisible=false, xticklabelsvisible=false, xlabelvisible=false,binsize=0.1)
        plot_raster_and_psth!(lg33, allcelldirs[cell_examples.both];previous=false, ylabelvisible=false, binsize=0.1)
        # tuning strength of individual cell vs contribution to decoder
        # plot performance on current poster vs previous poster
        lg1 = GridLayout(fig[2,1])
        lg11 = GridLayout(lg1[1,1])
        Label(lg11[1,1,TopLeft()],"D")
        ax = Axis(lg11[1,1])
        sc = scatter!(ax, μ_perf_c, μ_perf_p,color=ncells, markersize=30px)
        Colorbar(lg11[1,2], sc, label="No cells")
        ablines!(ax, 0.0, 1.0, linestyle=:dot, color=:black)
        ax.xlabel = "Performance current"
        ax.ylabel = "Performance previous"

        lg12 = GridLayout(lg1[1,2])
        ax41 = Axis(lg12[1,1])
        ax4 = Axis(lg12[2,1])
        linkxaxes!(ax4,ax41)
        ax41.xticklabelsvisible = false
        hist!(ax41, 1.0./current_poster_selectivity)
        scatter!(ax4, 1.0./current_poster_selectivity, dropdims(mean(contrib_current,dims=2),dims=2))
        ax4.xlabel = "Current poster\nselectivity"
        ax4.ylabel = "Relative coding contrib"
        vlines!(ax4, 1.0, linestyle=:dot, color=:black)
        vlines!(ax41, 1.0, linestyle=:dot, color=:black)
        rowsize!(lg12, 1, Relative(0.4))
        Label(lg12[1,1,TopLeft()], "E")

        ax51 = Axis(lg12[1,2])
        ax5 = Axis(lg12[2,2])
        Label(lg12[1,2,TopLeft()], "F")
        linkxaxes!(ax5,ax51)
        ax51.xticklabelsvisible = false
        hist!(ax51, 1.0./previous_poster_selectivity)
        scatter!(ax5, 1.0./previous_poster_selectivity, dropdims(mean(contrib_previous,dims=2),dims=2))
        vlines!(ax5, 1.0, linestyle=:dot, color=:black)
        vlines!(ax51, 1.0, linestyle=:dot, color=:black)
        ax5.xlabel = "Previous poster\nselectivity"
        rowsize!(fig.layout, 1, Relative(0.6))
        resize_to_layout!(fig)
        rowgap!(fig.layout, 1, 0)
        fig
    end
end

"""
Plot the view and place field for 5 cells with both place and view selectivity
"""
function plot_place_and_view_selective_cells()
    place_selective_cells = open("data/place_selective_cells.txt") do fid
       readlines(fid)
    end
    view_selective_cells = open("data/view_selective_cells.txt") do fid
       readlines(fid)
    end
    place_and_view_selective_cells = intersect(place_selective_cells, view_selective_cells)
    vm1,spm1 = cd(place_and_view_selective_cells[3]) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end
    vm2,spm2 = cd(place_and_view_selective_cells[4]) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end
     vm3,spm3 = cd(place_and_view_selective_cells[5]) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end
    vm4,spm4 = cd(place_and_view_selective_cells[6]) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end


    with_theme(poster_theme) do
        fig = Figure(size=(800,400))
        lg1 = GridLayout(fig[1,1])
        plot_view_and_place_fields!(lg1, vm1, spm1)
        lg2 = GridLayout(fig[1,2])
        plot_view_and_place_fields!(lg2, vm2, spm2;colorbar_label="")
        lg3 = GridLayout(fig[1,3])
        plot_view_and_place_fields!(lg3, vm3, spm3;colorbar_label="")
        lg4 = GridLayout(fig[1,4])
        plot_view_and_place_fields!(lg4, vm4, spm4;colorbar_label="")
        fig
    end
end

function plot_view_and_place_fields(celldir::String)
    vm,spm = cd(celldir) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end
    with_theme(poster_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_view_and_place_fields!(lg, vm, spm)
        fig
    end
end

function plot_view_and_place_fields!(lg, vm::ViewMapNew{T1}, spm::SpatialMapNew{T2};colorbar_label="Firing rate [Hz]") where T1 <: Real where T2 <: Real
    mm = get_maze_mesh()
    D = distancematrix(mm)
    m_floor = floor_topology3()
    D_floor = distancematrix(m_floor)
    f_sp = get_rate_map(spm;invalidate_unvisited=false)
    f_v = get_rate_map(vm;invalidate_unvisited=false)

    #smoothing
    Zp = fill_in_neighbours2(f_sp, D_floor, 12, T2(4.0))
    Zv = fill_in_neighbours2(vec(f_v), D, 12, T1(4.0))

    lscene = LScene(lg[1,1])
    plotmesh!(lscene, mm;floor_offset=-20, ceiling_offset=10, color=Zv,showsegments=false, colormap=:jet)
    m_floor2 = Translate(0.0, 0.0, -30)(m_floor)
    viz!(lscene, m_floor2;color=Zp,showsegments=false)
    lg2 = GridLayout(lg[1,2])
    Colorbar(lg2[1,1], colorrange=extrema(Zv), colormap=:jet, label=colorbar_label)
    Colorbar(lg2[2,1], colorrange=extrema(Zp), label=colorbar_label)
    lg
end

function plot_view_fields(patches::Vector{Vector{Vector{Int64}}})
    mm = Hippocampus.get_maze_mesh()
    Z,alpha = set_peaks(patches[1],mm,1.0)
    for (i,patch) in enumerate(patches[2:end])
        set_peaks!(Z, alpha, patch, i+!)
    end
    fig = explore(mm;color=Z, alpha=alpha,showsegments=true, floor_offset=-20, ceiling_offset=10.0,colormap=:brg)
end

function plot_decoding_results()
    #view coding
    data_view = JLD2.load(joinpath(@__DIR__, "..","data","view_decoding_preliminary_shuffled_results.jld2"))
    data_place = JLD2.load(joinpath(@__DIR__,"..","data","place_decoding_preliminary_results.jld2"))

    with_theme(poster_theme) do
        fig = Figure(size=(1448,698))
        lgp = GridLayout(fig[1,1])
        Label(lgp[1,1,TopLeft()],"A")
        plot_spatial_decoding_results!(lgp, data_place["km_results"], data_place["mean_err"], data_place["mean_err_sh"];_plot_theme=poster_theme)
        lgv = GridLayout(fig[1,2])
        Label(lgv[1,1,TopLeft()], "B")
        plot_view_decoding_results!(lgv, data_view["km_results"], data_view["mean_err"], data_view["mean_err_sh"];_plot_theme=poster_theme)
        fig
    end
end

function plot_knn_population_decoding_results(;kwargs...)
    fname = joinpath(@__DIR__, "..","data","categorical_decoding_run_more_cells.jld2")
    plot_knn_population_decoding_results(fname;kwargs...)
end

function plot_knn_population_decoding_results(fname::String;show_f1_score=false,_plot_theme=poster_theme,figsize=nothing)
    data = JLD2.load(fname)
    unique_categories = data["unique_categories"]
    mm = get_maze_mesh(;nrefinements=0)
    # categorise into floor, ceiling, walls and pillars
    tidx = categorize(mm)
    m_floor = floor_topology3(;nrefinements=0)
    m_floor2, m_ceiling, m_middle = Hippocampus.get_floor_and_ceiling(mm)
    # get the view categories
    Y = data["Y"]
    kidx_v = categorize(Y[1:3,:], mm) 
    # get the place cateogires
    kidx_p = mapto(Shadow("xy")(m_floor), Tuple.(eachcol(Y[4:5,:])))

    # floor has a centroid z-coordinate of 0, ceiling has a centroid z coordinate of 5
    # pillars have x,y centroid x,y coordinate larger than -12 and less than 12
    category = collect(zip(kidx_v, first.(kidx_p)))

    perf = data["perf"]
    if "fp_rate" in keys(data)
        fp_rate = data["fp_rate"]
    else
        fp_rate = zeros(size(perf)...)
    end
    if "fn_rate" in keys(data)
        fn_rate = data["fn_rate"]
    else
        fn_rate = zeros(size(perf)...)
    end
    f1_score = 2*perf./(2*perf .+ fp_rate .+ fn_rate)
    view_idx = [_cat[1] for _cat in unique_categories]
    place_idx = [_cat[2] for _cat in unique_categories]
    perf_view = zeros(maximum(view_idx))
    perf_place = zeros(maximum(place_idx))
    fp_view = fill!(similar(perf_view), 0.0)
    fn_view = fill!(similar(perf_view), 0.0)
    fp_place = fill!(similar(perf_place), 0.0)
    fn_place = fill!(similar(perf_place), 0.0)

    n_view = zeros(maximum(view_idx))
    n_place = zeros(maximum(place_idx))
    #TODO: Is this correct? Just taking the mean performance for each view, essentially
    # I think we need to take the priors into account as well
    prob_view = zeros(size(perf_view,1))
    prob_place= zeros(size(perf_place,1))
    for _cat in category
        (v,p) = _cat
        prob_view[v] += 1.0
        prob_place[p] += 1.0
    end
    prob_view ./= sum(prob_view)
    prob_place ./= sum(prob_place)
    for (j,_cat) in enumerate(unique_categories)
        vidx = _cat[1]
        pidx = _cat[2]
        fidx = isfinite.(perf[j,:])
        perf_view[vidx] += sum(perf[j,fidx]).*prob_place[pidx]
        fp_view[vidx] += sum(fp_rate[j,fidx]).*prob_place[pidx]
        fn_view[vidx] += sum(fn_rate[j,fidx]).*prob_place[pidx]
        perf_place[pidx] += sum(perf[j, fidx]).*prob_view[vidx]
        fp_place[pidx] += sum(fp_rate[j,fidx]).*prob_view[vidx]
        fn_place[pidx] += sum(fn_rate[j,fidx]).*prob_view[vidx]

        n_view[vidx] += sum(fidx).*prob_place[pidx]
        n_place[pidx] += sum(fidx).*prob_view[vidx]
    end
    f1_place = 2*perf_place./(2*perf_place + fp_place + fn_place)
    f1_view = 2*perf_view./(2*perf_view + fp_view + fn_view)
    perf_view ./= n_view
    perf_place ./= n_place

    # create view decoding conditioned on place
    pidx = sortperm(perf_place,rev=true)
    perf_view_place = zeros(size(perf_view,1), 6)
    f1_view_place = zeros(size(perf_view,1), 6)
    n_view_place = zeros(size(perf_view,1))
    fp_view_place = zeros(size(perf_view,1))
    fn_view_place = zeros(size(perf_view,1))
    for (i,_pidx) in enumerate(pidx[1:size(perf_view_place,2)])
        fill!(n_view_place, 0.0)
        fill!(fn_view_place, 0.0)
        fill!(fp_view_place, 0.0)
        _vidx = findall([_cat[2]==_pidx for _cat in unique_categories])
        view_idx = [_cat[1] for _cat in unique_categories[_vidx]]
        for (k,v) in enumerate(_vidx)
            fidx = isfinite.(perf[v,:])
            perf_view_place[view_idx[k],i] = sum(perf[v,fidx])
            fp_view_place[view_idx[k]] = sum(fp_rate[v,fidx])
            fn_view_place[view_idx[k]] = sum(fn_rate[v,fidx])
            n_view_place[view_idx[k]] = sum(fidx)
        end
        f1_view_place[:,i] = 2*perf_view_place[:,i]./(2*perf_view_place[:,i] .+ fp_view_place .+ fn_view_place)
        perf_view_place[:,i] ./= n_view_place
    end

    vidx = sortperm(perf_view, rev=true)
    perf_place_view = zeros(size(perf_place,1), 6)
    f1_place_view = zeros(size(perf_place,1), 6)
    n_place_view = zeros(size(perf_place,1))
    fp_place_view = zeros(size(perf_place,1))
    fn_place_view = zeros(size(perf_place,1))
    for (i,_vidx) in enumerate(vidx[1:size(perf_view_place,2)])
        fill!(n_place_view, 0.0)
        fill!(fn_place_view, 0.0)
        fill!(fp_place_view, 0.0)
        _pidx = findall([_cat[1]==_vidx for _cat in unique_categories])
        place_idx = [_cat[2] for _cat in unique_categories[_pidx]]
        for (k,v) in enumerate(_pidx)
            fidx = isfinite.(perf[v,:])
            perf_place_view[place_idx[k],i] = sum(perf[v,fidx])
            fp_place_view[place_idx[k]] = sum(fp_rate[v,fidx])
            fn_place_view[place_idx[k]] = sum(fn_rate[v,fidx])
            n_place_view[place_idx[k]] = sum(fidx)
        end
        f1_place_view[:,i] = 2*perf_place_view[:,i]./(2*perf_place_view[:,i] .+ fp_place_view .+ fn_place_view)
        perf_place_view[:,i] ./= n_place_view
    end

    #cluster in pca space
    X = data["X"]
    Xt,cat_t = generate_pseudosamples(X, category)
    # TODO: Place conditioned via and view-conditioned place
    xidx_p = findall(cc->cc[1].===vidx[1], cat_t)
    cat_p = [cc[2] for cc in cat_t[xidx_p]]
    xidx_v = findall(cc->cc[2].==pidx[1], cat_t)
    cat_v = [cc[1] for cc in cat_t[xidx_v]]

    pca = fit(PCA, Xt)
    Z = predict(pca, Xt)

    m_floor = Meshes.Translate(0.0, 0.0, -30)(m_floor)

    with_theme(_plot_theme) do
        if figsize === nothing
            figsize = (1596,1229)
        end
        fig = Figure(size=figsize)
        lg1 = GridLayout(fig[1,1])
        lgp1 = GridLayout(lg1[1,1])
        if show_f1_score
            _colorv = f1_view[tidx]
            _colorp = f1_place
            _label = "F1-score"
        else
            _colorv = perf_view[tidx]
            _colorp = perf_place
            _label = "Performance"
        end
        Label(lgp1[1,1,TopLeft()], "A")

        if plot_joint_matrix
            _zidx = CartesianIndex{2}.(unique_categories)
            ZZ = zeros(22,21)
            ZZ[_zidx] .= dropdims(mean(f1_score,dims=2),dims=2)
            axq = Axis(lgp1[1,1])
            hh = heatmap!(axq, ZZ, colormap=:Reds)
            Colorbar(lgp1[1,2], hh, label="F1-score")
            rowsize!(lg1, 1, Relative(0.4))
            axq.xlabel = "View bin"
            axq.ylabel = "Place bin"
        else
            lscene = LScene(lgp1[1,1],show_axis=false)
            plotmesh!(lscene, mm;color=_colorv, ceiling_offset=10, floor_offset=-10,colormap=:Purples,showsegments=true)
            viz!(lscene, m_floor;color=_colorp, colormap=:Greens, showsegments=true)
            lg12 = GridLayout(lgp1[1,2])
            Colorbar(lg12[1,1], colorrange=extrema(filter(isfinite,_colorv)), colormap=:Purples, label="$_label\nview",ticks=WilkinsonTicks(3))
            Colorbar(lg12[2,1], colorrange=extrema(filter(isfinite,_colorp)), colormap=:Greens, label="$_label\nplace",ticks=WilkinsonTicks(3))
        end
        # view probability conditioned on place
        if show_f1_score
            cr = extrema(filter(isfinite, f1_view_place))
        else
            cr = extrema(filter(isfinite, perf_view_place))
        end
        # cluster plot
        # place-conditioned view
        lgcc = GridLayout(lg1[2,1])
        Label(lgcc[1,1,TopLeft()],"B")
        # TODO: One for view-conditioned space, one for space conditioned view
        lscenep = LScene(lgcc[1,1],show_axis=false)
        _colors = HSV.(get_maze_category_colors())
        plotmesh!(lscenep, mm;color=_colors[tidx], ceiling_offset=10, floor_offset=-15, showsegments=true)
        floor_color = fill(parse(Colorant, :lightgray), nelements(m_floor))
        floor_color[pidx[1]] = parse(Colorant, :red)
        viz!(lscenep, m_floor;color=floor_color)
        lscenepc = Axis3(lgcc[2,1],xticklabelsvisible=false, yticklabelsvisible=false, zticklabelsvisible=false,
                                  xgridvisible=true, ygridvisible=true,zgridvisible=true,
                                  xlabelvisible=false, ylabelvisible=false, zlabelvisible=false, viewmode=:stretch)
        scatter!(lscenepc, Point3f.(eachcol(Z[1:3,xidx_v])),color=_colors[cat_v])

        # view-conditioned place
         Label(lgcc[1,2,TopLeft()],"C")
        # TODO: One for view-conditioned space, one for space conditioned view
        lscenev = LScene(lgcc[1,2],show_axis=false)
        _colors = fill(parse(Colorant, :lightgray), 22)
        _colors[vidx[1]]  = parse(Colorant, :red)
        plotmesh!(lscenev, mm;color=_colors[tidx], ceiling_offset=10, floor_offset=-15, showsegments=true)
        pcolor = resample_cmap(:tab20, 21)
        viz!(lscenev, m_floor;color=pcolor)
        lscenevc = Axis3(lgcc[2,2],xticklabelsvisible=false, yticklabelsvisible=false, zticklabelsvisible=false,
                                  xgridvisible=true, ygridvisible=true,zgridvisible=true,
                                  xlabelvisible=false, ylabelvisible=false, zlabelvisible=false, viewmode=:stretch)

        scatter!(lscenevc, Point3f.(eachcol(Z[1:3,xidx_p])),color=pcolor[cat_p])

        lg4 = GridLayout(fig[1,2])
        lscenes = [LScene(lg4[r,c], show_axis=false) for (r,c) in [(1,1),(1,2),(1,3),(2,1),(2,2),(2,3)]]
        Label(lg4[1,1,TopLeft()],"D")
        for k in 1:size(perf_view_place,2)
            lscene2 = lscenes[k]
            if show_f1_score
                _color = f1_view_place[:,k][tidx]
            else
                _color = perf_view_place[:,k][tidx]
            end
            plotmesh!(lscene2, mm;color=_color, ceiling_offset=10, floor_offset=-10,colormap=:Purples,colorrange=cr, showsegments=true)
            fcolor = fill(parse(Colorant, :lightgray), nelements(m_floor))
            fcolor[pidx[k]] = parse(Colorant, :red)
            viz!(lscene2, m_floor;color=fcolor, showsegments=true)
        end
        Colorbar(lg4[2,4], colorrange=cr, colormap=:Purples)

        # place probability conditioned on view
        if show_f1_score
            cr = extrema(filter(isfinite, f1_place_view))
        else
            cr = extrema(filter(isfinite, perf_place_view))
        end
        Label(lg4[3,1,TopLeft()],"E")
        lscenes = [LScene(lg4[r,c], show_axis=false) for (r,c) in [(3,1),(3,2),(3,3),(4,1),(4,2),(4,3)]]
        for k in 1:size(perf_place_view,2)
            lscene2 = lscenes[k]
            if show_f1_score
                _color = f1_place_view[:,k]
            else
                _color = perf_place_view[:,k]
            end
            mcolor = fill(parse(Colorant, :lightgray), nelements(mm))
            mcolor[tidx.==vidx[k]] .= parse(Colorant, :red)
            plotmesh!(lscene2, mm;color=mcolor, ceiling_offset=10, floor_offset=-10, showsegments=true)
            viz!(lscene2, m_floor;color=_color, showsegments=true,colorrange=cr, colormap=:Greens)
        end
        Colorbar(lg4[4,4], colorrange=cr, colormap=:Greens)
        colsize!(fig.layout, 1, Relative(0.3))
        fig
    end
end