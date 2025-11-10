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
        width = 53*72/2.5
        height = width
        fig = Figure(size=(width, height))
        # indivvidual cell responses
        lg2 = GridLayout(fig[1,1])

        lg21 = GridLayout(lg2[1,1])
        Label(lg21[1,1,TopLeft()], "A")
        plot_raster_and_psth!(lg21, allcelldirs[cell_examples.previous_poster_selective];previous=true,xlabelvisible=true, xticklabelsvisible=true,binsize=0.1,show_psth=show_psth)
        rowsize!(lg21, 1, Relative(0.6))
        lg31 = GridLayout(lg2[2,1])
        rowsize!(lg31, 1, Relative(0.6))
        plot_raster_and_psth!(lg31, allcelldirs[cell_examples.previous_poster_selective];previous=false, binsize=0.1, show_psth=show_psth)

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
        plot_raster_and_psth!(lg22, allcelldirs[cell_examples.poster_selective];previous=true,ylabelvisible=false,xlabelvisible=true, xticklabelsvisible=true,binsize=0.1, show_psth=show_psth)
        plot_raster_and_psth!(lg32, allcelldirs[cell_examples.poster_selective];previous=false, ylabelvisible=false, binsize=0.1, show_psth=show_psth)

        plot_raster_and_psth!(lg23, allcelldirs[cell_examples.both];previous=true, ylabelvisible=false, xticklabelsvisible=true, xlabelvisible=true,binsize=0.1, show_psth=show_psth)
        plot_raster_and_psth!(lg33, allcelldirs[cell_examples.both];previous=false, ylabelvisible=false, binsize=0.1,show_psth=show_psth)
        # tuning strength of individual cell vs contribution to decoder
        # plot performance on current poster vs previous poster
        lg1 = GridLayout(fig[2,1])
        colsize!(lg1, 1, Relative(0.4))
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
        ax41.ylabel = "# cells"
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
        # panel to show marginal histogram of coding cotribution
        lgpp = GridLayout(lg12[2,3])
        ax7 = Axis(lgpp[1,1], xticks=WilkinsonTicks(2))
        ax8 = Axis(lgpp[1,2], xticks=WilkinsonTicks(2))
        ax8.yticklabelsvisible = false
        hist!(ax7, dropdims(mean(contrib_current,dims=2),dims=2),direction=:x)
        hist!(ax8, dropdims(mean(contrib_previous,dims=2),dims=2), direction=:x)
        linkyaxes!(ax7, ax5, ax8)
        ax7.xlabel = "# cells"
        ax8.yticklabelsvisible = false
        rowsize!(fig.layout, 1, Relative(0.6))
        resize_to_layout!(fig)
        rowgap!(fig.layout, 1, 0)
        fig
    end
end

"""
Plot the view and place field for 5 cells with both place and view selectivity
"""
function plot_place_and_view_selective_cells(;_plot_theme=poster_theme)
    place_selective_cells = open("data/place_selective_cells.txt") do fid
       readlines(fid)
    end
    view_selective_cells = open("data/view_selective_cells.txt") do fid
       readlines(fid)
    end
    place_and_view_selective_cells = intersect(place_selective_cells, view_selective_cells)
    place_or_view_selective_cells = union(place_selective_cells, view_selective_cells)
    svm_spm = map(place_or_view_selective_cells[[4,13,21,9]]) do celldir
        cd(celldir) do
            svm = Hippocampus.SmoothedMap(Hippocampus.ViewMapNew, Hippocampus.UnityRaytraceData;min_speed=2.0, trial_start=2,method=:adaptive, α=1000.0^2, rmax=10)
            spm = Hippocampus.SmoothedMap(Hippocampus.SpatialMapNew;min_speed=2.0, trial_start=2,method=:adaptive, α=10000.0^2, rmax=10) 
            svm,spm
        end
    end

    with_theme(_plot_theme) do
        fig = Figure(size=(800,400))
        lg1 = GridLayout(fig[1,1])
        plot_view_and_place_fields!(lg1, svm_spm[1]...)
        lg2 = GridLayout(fig[1,2])
        plot_view_and_place_fields!(lg2, svm_spm[2]...;colorbar_label="")
        lg3 = GridLayout(fig[1,3])
        plot_view_and_place_fields!(lg3, svm_spm[3]...;colorbar_label="")
        lg4 = GridLayout(fig[1,4])
        plot_view_and_place_fields!(lg4, svm_spm[4]...;colorbar_label="")
        fig
    end
end

function plot_view_and_place_fields(celldir::String;kwargs...)
    vm,spm = cd(celldir) do
        vm = Hippocampus.ViewMapNew(Hippocampus.UnityRaytraceData;min_speed=2.0,trial_start=2,do_save=false)
        spm = Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
        vm,spm
    end
    svm, spm = cd(celldir) do
        svm = Hippocampus.SmoothedMap(Hippocampus.ViewMapNew, Hippocampus.UnityRaytraceData;min_speed=2.0, trial_start=2,method=:adaptive, α=1000.0^2, rmax=10)
        spm = Hippocampus.SmoothedMap(Hippocampus.SpatialMapNew;min_speed=2.0, trial_start=2,method=:adaptive, α=10000.0^2, rmax=10)
        svm, spm
    end
    with_theme(poster_theme) do
        fig = Figure(size=(700, 700))
        lg = GridLayout(fig[1,1])
        plot_view_and_place_fields!(lg, svm, spm;kwargs...)
        fig
    end
end

function plot_view_and_place_fields!(lg, vm::ViewMapNew{T1}, spm::SpatialMapNew{T2};σp=T2(4.0), σv=T1(σp), colorbar_label="Firing rate [Hz]") where T1 <: Real where T2 <: Real
    mm = get_maze_mesh()
    D = distancematrix(mm)
    m_floor = floor_topology3()
    m_floor_flat = Shadow("xy")(m_floor)
    D_floor = distancematrix(m_floor_flat)
    f_sp = get_rate_map(spm;invalidate_unvisited=false)
    f_v = get_rate_map(vm;invalidate_unvisited=false)

    #smoothing
    #Zp = fill_in_neighbours2(f_sp, D_floor, round(Int64, 3*σp), σp)
    #Zv = fill_in_neighbours2(vec(f_v), D, round(Int64, 3*σv), σv)
    Z_spm_s, X_spm_s, Y_spm_s = adaptive_smoothing(spm.weight, spm.occupancy, m_floor_flat, (10000.0)^2,rmax=10)
    Y_spm_s[spm.occupancy.==0.0] .= 0.0
    Z_spm_s[spm.occupancy .==0.0] .= NaN
    Z_vm_s, X_vm_s, Y_vm_s = adaptive_smoothing(vm.weight, vm.occupancy, mm, (1000.0)^2,rmax=10)
    Y_vm_s[vm.occupancy.==0.0] .= 0.0
    #Z_vm_s[vm.occupancy .==0.0] .= NaN


    sic_place = compute_sic(X_spm_s, Y_spm_s)
    sic_view = compute_sic(X_vm_s, Y_vm_s)

    lscene = LScene(lg[1,1])
    plotmesh!(lscene, mm;floor_offset=-15, ceiling_offset=10, color=Z_vm_s,showsegments=false, colormap=:jet)
    m_floor2 = Translate(0.0, 0.0, -30)(m_floor)
    viz!(lscene, m_floor2;color=Z_spm_s,showsegments=false)
    # labels
    text!(lscene, 0.0, 0.9, text="SIC = $(round(sic_view; sigdigits=2))", space=:relative,
                            rotation=-π/2)
    text!(lscene, 0.0, 0.3, text="SIC = $(round(sic_place;sigdigits=2))", space=:relative,
                            rotation=-π/2)
    lg2 = GridLayout(lg[1,2])
    Colorbar(lg2[1,1], colorrange=extrema(filter(isfinite,Z_vm_s)), colormap=:jet, label=colorbar_label)
    Colorbar(lg2[2,1], colorrange=extrema(filter(isfinite, Z_spm_s)), label=colorbar_label)
    lg
end

function plot_view_and_place_fields!(lg, svm::SmoothedMap, sspm::SmoothedMap;colorbar_label="Firing rate [Hz]")
    mm = get_maze_mesh()
    mm_simple = get_maze_mesh(;nrefinements=0)
    m_floor = floor_topology3()
    m_floor_simple = floor_topology3(;nrefinements=0)
    m_floor_flat = Shadow("xy")(m_floor)

    sic_place = compute_skaggs_sic(sspm)
    sic_view = compute_skaggs_sic(svm)
    Z_vm_s = get_rate_map(svm;invalidate_unvisited=false)
    Z_spm_s = get_rate_map(sspm;invalidate_unvisited=true)

    lscene = LScene(lg[1,2],show_axis=false)
    plotmesh!(lscene, mm;floor_offset=-15, ceiling_offset=10, color=Z_vm_s,showsegments=false, colormap=:jet)
    plotmesh!(lscene, mm_simple;floor_offset=-15, ceiling_offset=10, alpha=0,showsegments=true, colormap=:jet)
    m_floor2 = Translate(0.0, 0.0, -40)(m_floor)
    viz!(lscene, m_floor2;color=Z_spm_s,showsegments=false)
    viz!(lscene, Translate(0.0, 0.0, -40)(m_floor_simple);alpha=0,showsegments=true)
    # indicate north
    arrows3d!(lscene, Point3f(0.0, 10.0, -40.0), Point3f(0.0, 3.0, 0.0),color=:black)
    # labels
    text!(lscene, 0.9, 0.9, text="SIC = $(round(sic_view; sigdigits=2))", space=:relative,
                            rotation=-π/2)
    text!(lscene, 0.9, 0.3, text="SIC = $(round(sic_place;sigdigits=2))", space=:relative,
                            rotation=-π/2)
    lg2 = GridLayout(lg[1,1])
    Colorbar(lg2[1,1], colorrange=extrema(filter(isfinite,Z_vm_s)), colormap=:jet, label=colorbar_label,flipaxis=false)
    Colorbar(lg2[2,1], colorrange=extrema(filter(isfinite, Z_spm_s)), label=colorbar_label,flipaxis=false)
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
        #lgp1 = GridLayout(lg1[2,1])
        if show_f1_score
            _colorv = f1_view[tidx]
            _colorp = f1_place
            _label = "F1-score"
        else
            _colorv = perf_view[tidx]
            _colorp = perf_place
            _label = "Performance"
        end
        Label(lg1[2,1,TopLeft()], "C")

        if !isempty(fname_ind)
           lgpq = GridLayout(lg1[2,1])
           plot_independent_vs_joint_category_decoding!(lgpq, fname_ind, fname;_plot_theme=_plot_theme) 
        elseif plot_joint_matrix
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
        lgcc = GridLayout(lg1[1,1])
        Label(lgcc[1,1,TopLeft()],"A")
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
         Label(lgcc[1,2,TopLeft()],"B")
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


function plot_independent_vs_joint_category_decoding()
    plot_independent_vs_joint_category_decoding("data/categorical_decoding_run_more_cells_pca_f1_score_train_independent_test_joint.jld2", "data/categorical_decoding_run_more_cells_pca_f1_score_train_joint_test_joint.jld2")
end

function plot_independent_vs_joint_category_decoding(fname_independent::String, fname_joint::String;_plot_theme=poster_theme,kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(640,554))
        lg = GridLayout(fig[1,1])
        plot_independent_vs_joint_category_decoding!(lg, fname_independent, fname_joint;_plot_theme=_plot_theme,kwargs...)
        fig
    end
end

function plot_independent_vs_joint_category_decoding!(lg, fname_independent::String, fname_joint::String;_plot_theme=poster_theme,nshuffles=1000)
    ind_data = JLD2.load(fname_independent)
    joint_data = JLD2.load(fname_joint)

    # compute f1-score for bot
    f1_score_ind = 2*ind_data["perf"]./(2*ind_data["perf"] .+ ind_data["fp_rate"] .+ ind_data["fn_rate"])
    f1_score_joint= 2*joint_data["perf"]./(2*joint_data["perf"] .+ joint_data["fp_rate"] .+ joint_data["fn_rate"])

    f1_score_ind_mean = zeros(size(f1_score_ind,1))
    for i in axes(f1_score_ind,1)
        fidx = isfinite.(f1_score_ind[i,:])
        f1_score_ind_mean[i] = mean(f1_score_ind[i,fidx])
    end

    f1_score_joint_mean = zeros(size(f1_score_joint,1))
    for i in axes(f1_score_joint,1)
        fidx = isfinite.(f1_score_joint[i,:])
        f1_score_joint_mean[i] = mean(f1_score_joint[i,fidx])
    end
    
    #make sure we are using the same categories
    unique_cat_ind = ind_data["unique_categories"]
    unique_cat_joint = joint_data["unique_categories"]
    unique_categories = intersect(unique_cat_ind, unique_cat_joint)
    qidx_ind = findall(in(unique_categories).(unique_cat_ind))
    qidx_joint = findall(in(unique_categories).(unique_cat_joint))

    #make sure they match up
    qidx2_ind = [findfirst(qq->cc==qq, unique_categories) for cc in unique_cat_ind[qidx_ind]]
    qidx2_joint = [findfirst(qq->cc==qq,unique_categories) for cc in unique_cat_joint[qidx_joint]]
    qidx_ind = qidx_ind[qidx2_ind]
    qidx_joint = qidx_joint[qidx2_joint]

    f1_score_joint_mean = f1_score_joint_mean[qidx_joint]
    f1_score_ind_mean = f1_score_ind_mean[qidx_ind]
    #shuffle test
    nq = sum(f1_score_joint_mean .> f1_score_ind_mean)
    nqs = zeros(nshuffles)
    f1_score_joint_mean_sh = fill!(similar(f1_score_joint_mean), 0.0)
    f1_score_ind_mean_sh = fill!(similar(f1_score_ind_mean), 0.0)
    f1_score_mean = [f1_score_joint_mean f1_score_ind_mean]
    for i in 1:nshuffles
        for j in 1:length(f1_score_joint_mean_sh)
            i1,i2 = shuffle(1:2)
            f1_score_joint_mean_sh[j] = f1_score_mean[j,i1]
            f1_score_ind_mean_sh[j] = f1_score_mean[j,i2]
        end
        nqs[i] = sum(f1_score_joint_mean_sh .> f1_score_ind_mean_sh)
    end
    @show nq percentile(nqs, 99)
    with_theme(_plot_theme) do
        ax = Axis(lg[1,1])
        scatter!(ax,f1_score_joint_mean, f1_score_ind_mean)
        ablines!(ax, [0.0], [1.0], linestyle=:dot, color=:black)
        ax.xticks = WilkinsonTicks(3)
        ax.yticks = WilkinsonTicks(3)
        ax.xlabel = "F1-score joint"
        ax.ylabel = "F1-score independent"
    end
end

function plot_view_and_place_decoding(fname_place::String, fnane_view::String;_plot_theme=poster_theme)
    data_v = JLD2.load(fnane_view)
    data_p = JLD2.load(fname_place)

    f1_score_p = 2*data_p["perf"]./(2*data_p["perf"] .+ data_p["fp_rate"] .+ data_p["fn_rate"])
    f1_score_v = 2*data_v["perf"]./(2*data_v["perf"] .+ data_v["fp_rate"] .+ data_v["fn_rate"])

    mm = get_maze_mesh(;nrefinements=0)
    tidx = categorize(mm)
    m_floor = floor_topology3(;nrefinements=0)

    cat_view = data_v["unique_categories"]
    sidx = sortperm(cat_view)
    color_v = dropdims(mean(f1_score_v,dims=2),dims=2)[sidx][tidx]
    cat_place = data_p["unique_categories"]
    sidx = sortperm(cat_place)
    color_p = dropdims(mean(f1_score_p,dims=2),dims=2)[sidx]
    with_theme(_plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        lscene = LScene(lg[1,1], show_axis=false)
        plotmesh!(lscene, mm;color=color_v, ceiling_offset=10, floor_offset=-10,colormap=:Purples,showsegments=true)
        viz!(lscene, Translate(0.0, 0.0, -35.0)(m_floor);color=color_p, colormap=:Greens, showsegments=true)
        lg12 = GridLayout(lg[1,2])
        Colorbar(lg12[1,1], colorrange=extrema(filter(isfinite,color_v)), colormap=:Purples, label="F1-score\nview",ticks=WilkinsonTicks(3))
        Colorbar(lg12[2,1], colorrange=extrema(filter(isfinite,color_p)), colormap=:Greens, label="F1-scorel\nplace",ticks=WilkinsonTicks(3))
        fig
    end


end

function plot_occupancy(Y::Matrix{T},twin::Vector{T};_plot_theme=poster_theme) where T <: Real
    mm = get_maze_mesh(;nrefinements=0)

    m_floor = floor_topology3(;nrefinements=0)
    m_floor_flat = Shadow("xy")(m_floor)
    kidx_v = mapto(mm, Tuple.(eachcol(Y[1:3,:])))
    kidx_p = mapto(m_floor_flat, Tuple.(eachcol(Y[4:5,:])))
    duration_p = zeros(T, nelements(m_floor))
    V
    duration_v = zeros(T, maximum(tidx))
    num_v = zeros(T, maximum(tidx))
    num_p = zeros(T, nelements(m_floor)) 
    for (i,t) in enumerate(twin)
        duration_v[tidx[kidx_v[i]]] .+= t
        num_v[tidx[kidx_v[i]]] .+= 1.0 
        duration_p[kidx_p[i]] .+= t
        num_p[kidx_p[i]] .+= 1.0
    end
    duration_v ./= counts(tidx) 
    with_theme(_plot_theme) do
        fig = Figure(size=(700,700))
        lg = GridLayout(fig[1,1])
        lscene = LScene(lg[1,1])
        plotmesh!(lscene, mm;color=duration_v[tidx], ceiling_offset=10, floor_offset=-15, colormap=:Purples, showsegments=true)
        viz!(lscene, Translate(0.0, 0.0, -30)(m_floor);color=duration_p, colormap=:Greens)
        lgm = GridLayout(lg[1,2])
        Colorbar(lgm[1,1], colorrange=extrema(duration_v), colormap=:Purples, label="Mean duration\nview [s]")
        Colorbar(lgm[2,1], colorrange=extrema(duration_p), colormap=:Greens, label="Mean duration\nplace [s]")
        fig
    end
end


function plot_maze_with_posters(;kwargs...)
    with_theme(poster_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_maze_with_posters!(lg;kwargs...)
        fig
    end
end

function illustrate_allocentric_vs_egocentric()
    with_theme(poster_theme) do
        width = (40/2.5)*72
        height = 0.8*width
        fig = Figure(size=(width,height))
        lg1 = GridLayout(fig[1,1])
        Label(lg1[1,1,TopLeft()],"A")
        # 3D view of the maze
        plot_maze_with_posters!(lg1)
        lg2 = GridLayout(fig[1,2])
        #lg21 = GridLayout(lg2[1,1])
        # allocentric
        plot_flat_maze_with_posters!(lg2;start_point=Point2f(-5.0, -10.0), end_point=Point2f(6.0, 10.0))
        Label(lg2[1,1,TopLeft()], "B")
        #egocentric
        lg22 = GridLayout(fig[2,1:2])
        illustrate_viewpoint_coding!(lg22)
        Label(lg22[1,1,TopLeft()], "C")
        rowsize!(fig.layout, 1, Relative(0.4))
        fig
    end
end
function illustrate_viewpoint_coding()
    with_theme(poster_theme) do
        fig = Figure(size=(800, 600))
        lg = GridLayout(fig[1,1])
        illustrate_viewpoint_coding!(lg)
        fig
    end
end
"""
Illustrate how viewpoint coding gets you from the cat poster to the donkey poster
"""
function illustrate_viewpoint_coding!(lg)
    with_theme(poster_theme) do
        lg1 = GridLayout(lg[1,1])
        Label(lg1[1,1,Top()],"1")
        plot_maze_with_posters!(lg1;eyepos=Makie.Vec3f(-5.0, -10.0, 1.0), lookat=Makie.Vec3f(-4.5, -7.5, 1.0))
        lg2 = GridLayout(lg[1,2])
        Label(lg2[1,1,Top()],"2")
        plot_maze_with_posters!(lg2;eyepos=Makie.Vec3f(-3.0, -10.0, 1.0), lookat=Makie.Vec3f(5.0, 2.5, 1.0))
        lg3 = GridLayout(lg[2,1])
        Label(lg3[1,1,Top()],"3")
        plot_maze_with_posters!(lg3;eyepos=Makie.Vec3f(-3.0, 10.0, 1.0), lookat=Makie.Vec3f(5.0, 7.5, 1.0))
        lg4 = GridLayout(lg[2,2])
        Label(lg4[1,1,Top()],"4")
        plot_maze_with_posters!(lg4;eyepos=Makie.Vec3f(6.0, 10.0, 1.0), lookat=Makie.Vec3f(4.5, 7.5, 1.0))
    end
end

function plot_maze_with_posters!(lg;eyepos::Union{Makie.Vec3f,Nothing}=nothing, lookat::Union{Makie.Vec3f,Nothing}=nothing)
    mm_simple = get_maze_mesh(nrefinements=0)
    tidx = categorize(mm_simple)
    mm = get_maze_mesh()
    posters = Posters(mm, Hippocampus.poster_pos) 
    # fill in colors
    colors = fill(parse(Colorant, :lightgray),nelements(mm_simple))  
    colors[tidx.==2] .= HSV(parse(Colorant, :bisque3))
    colors[in([3,4,5,6]).(tidx)] .= HSV(parse(Colorant, :bisque2)) 
    colors[in([7,8,9,10]).(tidx)] .= HSV(212,70,100)
    colors[in([11,12,13,14]).(tidx)] .= HSV(119,70,100)
    colors[in([15,16,17,18]).(tidx)] .= HSV(356,70,100) 
    colors[in([19,20,21,22]).(tidx)] .= HSV(46,70,100)
    if eyepos !== nothing || lookat !== nothing
        hide_ceiling = false
    else
        hide_ceiling = true
    end
    with_theme(poster_theme) do
        lscene = LScene(lg[1,1],show_axis=false)  
        plot!(lscene, posters, shading=false)
        Hippocampus.plotmesh!(lscene, mm_simple;color=colors,showsegments=true, hide_ceiling=hide_ceiling, floor_offset=0.0)
        if eyepos !== nothing || lookat !== nothing
            cc = cam3d!(lscene.scene, center=false, fov=60) 
            update_cam!(lscene.scene, cc, eyepos, lookat)
        end
    end
end

function plot_flat_maze_with_posters()
    with_theme(poster_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_flat_maze_with_posters!(lg)
        fig
    end
end

function plot_flat_maze_with_posters!(lg;start_point::Union{Nothing, Point2f}=nothing, end_point::Union{Nothing, Point2f}=nothing)
    images = Dict(k=>load(v) for (k,v) in poster_img)
    colors = [HSV(46,70,100), #yellow
              HSV(212,70,100), #blue
              HSV(356,70,100), #red
              HSV(119,70,100) #geen
            ] 

    with_theme(poster_theme) do
        ax = Axis(lg[1,1],backgroundcolor=:bisque3)
        ax.xticksvisible = false
        ax.xticklabelsvisible = false
        ax.yticksvisible = false
        ax.yticklabelsvisible = false
        ax.topspinevisible = true
        ax.rightspinevisible = true

        limits!(ax, -12.5, 12.5, -12.5, 12.5)
        scatter!(ax, [Point2f(-5,5), Point2f(-5,-5), Point2f(5,5), Point2f(5,-5)], marker=Rect, markerspace=:data,color=colors, markersize=5)
        scatter!([Point2f(v[1:2]...) for (k,v) in poster_pos],  marker=[images[k] for (k,v) in poster_pos], markersize=4, markerspace=:data)
        if start_point !== nothing
            scatter!(ax, start_point, color=:gray)
        end
        if end_point !== nothing
            scatter!(ax, end_point, color=:black)
        end
        arrows!(ax, Point2f(0.0, 10.0), Point2f(0.0, 2.0), color=:black,arrowsize=10.0)
    end
end

function plot_gaze_path(unity_data::UnityRaytraceData,trialidx;_plot_theme=poster_theme)
    with_theme(_plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_gaze_path!(lg, unity_data, trialidx;_plot_theme=_plot_theme)
        fig
    end
end

function plot_gaze_path!(lg, unity_data::UnityRaytraceData, trialidx::Observable{Vector{Int64}}=Observable([1]);_plot_theme=poster_theme)
    mm = get_maze_mesh()
    D = distancematrix(mm)
    mm_simple = get_maze_mesh(;nrefinements=0)
    nt = numtrials(unity_data)
    points = lift(trialidx) do _ti
        traj = Point3f[]
        for ti in _ti 
            if 0 < ti <= nt 
                tt,tg,tp,fixmask,fo = get_trial(unity_data, ti;trial_start=2)
                # figure out which points are not on the maze walls
                kidx = Hippocampus.mapto(mm, Tuple.(eachcol(tg))) 
                fidx = findall((!isempty).(kidx))
                # we also want to split up the path, such that we do not connect across regions that not connected 
                qidx = Vector{Int64}[] 
                push!(qidx, [fidx[1]])
                for f in fidx[2:end] 
                    _kidx0 = first(kidx[qidx[end][end]])
                    _kidx1 = first(kidx[f])
                    if D[_kidx0, _kidx1] <= 2 
                        push!(qidx[end], f)
                    else
                        push!(qidx, [f])
                    end
                end
                for _qidx in qidx
                    append!(traj, Point3f.(eachcol(tg[:,_qidx])))
                    push!(traj, Point3f(NaN))
                end
            else
                traj = [Point3f(NaN)]
            end
        end
        traj
    end
    with_theme(_plot_theme) do
        lscene = LScene(lg[1,1],show_axis=false)
        plotmesh!(lscene, mm_simple;hide_ceiling=true, showsegments=true,alpha=0.0)
        lines!(lscene, points)
    end
end

"""
Group trials by pairs of poster ids
"""
function group_trials(triggers::Matrix{T}) where T<:Union{T2, Missing} where T2 <: Integer
    nt = size(triggers,1)
    groups = Dict{Tuple{Int64, Int64}, Vector{Int64}}()

    prev_posterid = 0
    for i in 1:nt
        # make sure the trial was correct
        if ismissing(triggers[i,3])
            continue
        end
        if !(30 < triggers[i,3] < 40)
            posterid = 0
        else
            posterid = triggers[i,1] - 10
            kk = (prev_posterid, posterid)
            if !(kk in keys(groups))
                groups[kk] = Int64[]
            end
            push!(groups[kk], i)
        end
        prev_posterid = posterid
    end
    groups
end

function run_lesion_simulation(;nruns=100,nbatch=20)
    # get all cells
    allcelldirs = open("/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/cell_list.txt") do fid
        readlines(fid)
    end
    vpr = map(allcelldirs) do celldir
             cd(celldir) do
                  try
                    return Hippocampus.ViewRepresentation(Hippocampus.UnityRaytraceData)
                  catch ee
                   return nothing
                   end

              end
    end;
    fidx = findall((!isnothing).(vpr))
    place_selective_cells = open("data/place_selective_cells.txt") do fid
        readlines(fid)
    end
    view_selective_cells = open("data/view_selective_cells.txt") do fid
       readlines(fid)
    end

    place_or_view_selective_cells = union(place_selective_cells, view_selective_cells)

    vidx = fidx[findall(in(place_or_view_selective_cells).(allcelldirs[fidx]))]
    nidx = fidx[findall((!in(place_or_view_selective_cells)).(allcelldirs[fidx]))]

    # 1) Remove all place or view selective cells
    no_place_or_view_selective_cells = allcelldirs[nidx]

    # 2) Remove an equivalent number that is not view or place selective
    f1_score_mean = Vector{Vector{Float64}}(undef, nbatch)
    unique_categories = Vector{Vector{Tuple{Int64, Int64}}}(undef, nbatch)
    for i in 1:nbatch
        pidx = setdiff(fidx, shuffle(nidx)[1:length(place_or_view_selective_cells)])
        with_place_and_view_selective_cells = allcelldirs[pidx]
        
        @assert length(no_place_or_view_selective_cells) == length(with_place_and_view_selective_cells)
        @assert isempty(intersect(no_place_or_view_selective_cells, place_or_view_selective_cells))
        @assert sort(intersect(with_place_and_view_selective_cells,place_or_view_selective_cells)) == sort(place_or_view_selective_cells)

        X,Y,twin = Hippocampus.get_population_representation(something.(vpr[pidx]))
        perf,fp_rate,fn_rate,unique_categories[i] = Hippocampus.population_decoder_simple(X,Y;k=30,nruns=nruns,do_pca=true,decode_view=true, decode_place=true,joint=true)
        f1_score = 2*perf./(2*perf .+ fp_rate .+ fn_rate)
        _f1_score_mean = zeros(size(f1_score,1))
        for i in axes(f1_score,1)
            _fidx = isfinite.(f1_score[i,:])
            _f1_score_mean[i] = mean(f1_score[i,_fidx])
        end
        f1_score_mean[i] = _f1_score_mean
    end
    f1_score_mean, unique_categories
end

function plot_lesion_results()
    fname = "data/categorical_decoding_run_remove_non_place_and_view_train_joint_test_joint_new.jld2"
    qq = JLD2.load(fname)
    fname2 = "data/categorical_decoding_run_no_place_or_view_cells_pca_f1_score_train_joint_test_joint.jld2"
    qq2 = JLD2.load(fname2)
    f1_score = 2*qq2["perf"]./(2*qq2["perf"] .+ qq2["fp_rate"] .+ qq2["fn_rate"])
    _f1_score_mean = zeros(size(f1_score,1))
    for i in axes(f1_score,1)
        _fidx = isfinite.(f1_score[i,:])
        _f1_score_mean[i] = mean(f1_score[i,_fidx])
    end

    f1_and_view = Dict()
    for (i,k) in enumerate(qq2["unique_categories"])
        a = f1_score[i]
        b = Float64[]
        for (j,(_f1_score, uq)) in enumerate(zip(qq["f1_scor_mean"], qq["unique_categories"]))
            idx = findfirst(cc->cc==k, uq)
            if idx !== nothing
                push!(b, _f1_score[idx])
            end
        end
        f1_and_view[k] = (a, b)
    end
    
    with_theme(poster_theme) do
        fig = Figure()
        ax = Axis(fig[1,1])
        xx = cat([f1_and_view[k][2] for k in keys(f1_and_view)]...,dims=1)
        yy = cat([fill(f1_and_view[k][1], length(f1_and_view[k][2])) for k in keys(f1_and_view)]...,dims=1)
        scatter!(ax, xx, yy)
        @show sum(yy.>xx)./length(xx)
        ablines!(ax, 0.0, 1.0, linestyle=:dot, color=:black)
        ax.xlabel = "F1 Remove non-place/view cells"
        ax.ylabel = "F1 Remove place-view cells"

        ax2 = Axis(fig[2,1])
        xx2 = 1:length(qq2["unique_categories"]) 
        yy3_l = Float64[]
        yy3_u = Float64[]
        xx3 = Float64[]
        @show qq2["unique_categories"]
        for k in keys(f1_and_view)
            _f1 = f1_and_view[k][2]
            u,l = percentile(_f1, [25,75])
            push!(yy3_u, u)
            push!(yy3_l, l)
            ii = findfirst(cc->cc=k, qq2["unique_categories"])
            push!(xx3, ii)
        end
        #yy3_l = cat([percentile(f1_and_view[k][2],25) for k in keys(f1_and_view)]...,dims=1)
        #yy3_u = cat([percentile(f1_and_view[k][2],75) for k in keys(f1_and_view)]...,dims=1)
        #xx3 = cat([fill(findfirst(cc->cc==k, qq2["unique_categories"]),length(f1_and_view[k][2])) for k in keys(f1_and_view)]...,dims=1)
        #xx3 = cat([findfirst(cc->cc==k, qq2["unique_categories"]) for k in keys(f1_and_view)]...,dims=1)
        scatter!(ax2, xx3, yy3_l, color=Cycled(1))
        scatter!(ax2, xx3, yy3_u,color=Cycled(1))
        scatter!(ax2, xx2, _f1_score_mean)
        fig
    end
end