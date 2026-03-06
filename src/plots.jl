using CairoMakie


function plot_place_map(celldir::String;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(1000,250))
        lg = GridLayout(fig[1,1])
        plot_place_map!(lg, celldir;kwargs...)
        fig
    end
end

"""
Illustrate place fields for the cell in `celldir`.

One panel with the raw firing rate map, one with a gaussian smoothed firing rate map and one with the adaptively smoothed firing rate map
"""
function plot_map!(lg, ::Type{T}, celldir::String;nrefinements=(p=3,g=2),σ=3, α=1000.0.^2, ylabel="", smooth_method::Vector{Symbol}=[:gaussian,:adaptive], smooth_params=[(σ=3,), (α=1000.0^2,)], kwargs...) where T <: AbstractInformationContent
    sm, smg, sma,mm, sic,sicg,sica = cd(celldir) do
        #sp = Spiketrain()
        #rp = cd(DPHT.process_level(level(RippleData))) do
        #    RippleData()
        #end
        #unity_gaze_data = cd(DPHT.process_level(level(UnityRaytraceData))) do
        #    UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
        #end
        #jocc = cd(DPHT.process_level(level(JointOccupancy))) do
        #    JointOccupancy(;redo=false, nrefinements=nrefinements,trial_start=trial_start)
        #end
        #jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
        sic = compute_skaggs_sic(T,10_000;kwargs...)
        mm = get_mesh(T,nrefinements)
        #mm = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        #vpvrpb = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...)
        jmb = JointMap(;kwargs...)
        # TODO: I need to be able to load this directly
        spmb = maptype(T)(jmb,mm)
        smg = SmoothedMap(spmb;method=smooth_method[1], smooth_params[1]...)
        sicg = compute_skaggs_sic(T,10_000;smooth=true, smoothing_method=smooth_method[1], smooth_params[1]...,kwargs...)
        sma = SmoothedMap(spmb;method=smooth_method[2], smooth_params[2]...)
        sica = compute_skaggs_sic(T,10_000;smooth=true, smoothing_method=smooth_method[2], smooth_params[2]...,kwargs...)
        spmb,smg,sma,mm, sic,sicg,sica
    end
    @show sic.sic0 percentile(sic.sic, 95)
    @show sicg.sic0 percentile(sicg.sic, 95)
    @show sica.sic0 percentile(sica.sic, 95)
    # add SIC with distribution
    ax0 = Axis(lg[1,1])
    boxplot!(ax0, fill(1.0, length(sic.sic)), sic.sic;show_outliers=false, show_notch=true,color=:gray)
    hlines!(ax0, percentile(sic.sic,95), linestyle=:dot)
    ax0.xticklabelsvisible = false
    ax0.xticksvisible = false
    ax0.bottomspinevisible = false
    scatter!(ax0, [1.0],[sic.sic0],color=:red)
    ax0.ylabel = "SIC"

    if embeddim(mm) == 2
        ax1 = Axis(lg[1,2])
        ax2 = Axis(lg[1,4])
        ax3 = Axis(lg[1,6])
    else
        ax1 = LScene(lg[1,2],show_axis=false)
        ax2 = LScene(lg[1,4], show_axis=false)
        ax3 = LScene(lg[1,6], show_axis=false)
    end

    Z = sm.weight./sm.occupancy
    if embeddim(mm) == 2
        viz!(ax1, mm;showsegments=false, color=:lightgray)
        viz!(ax1, mm;showsegments=false, color=Z)
        ax1.title = "Raw"
    else
        plotmesh!(ax1, mm;color=:lightgray, showsegments=false, floor_offset=-20, ceiling_offset=10)
        plotmesh!(ax1, mm;color=Z, showsegments=false, floor_offset=-20, ceiling_offset=10)
    end
    Colorbar(lg[1,3], colorrange=(extrema(filter(isfinite, Z))),label="Firing rate [Hz]")
    ax1.title = "Raw"
    Zg = smg.weight./smg.occupancy
    Zg[smg.unvisited] .= NaN 
    if embeddim(mm) == 2
        viz!(ax2, mm;showsegments=false, color=:lightgray)
        viz!(ax2, mm;showsegments=false, color=Zg)
        ax2.title =  strip(string(smooth_params[1]),['(',')'])
    else
        plotmesh!(ax2, mm;color=:lightgray, showsegments=false, floor_offset=-20, ceiling_offset=10)
        plotmesh!(ax2, mm;color=Zg, showsegments=false, floor_offset=-20, ceiling_offset=10)
    end
    Colorbar(lg[1,5], colorrange=(extrema(filter(isfinite, Zg))), label="Firing rate [Hz]")

    ax2.title = "σ = $(σ)"
    Za = sma.weight./sma.occupancy
    Za[sma.unvisited] .= NaN
    if embeddim(mm) == 2
        viz!(ax3, mm;showsegments=false, color=:lightgray)
        viz!(ax3, mm;showsegments=false, color=Za)
        ax3.title =  strip(string(smooth_params[2]),['(',')'])
    else
        plotmesh!(ax3, mm;color=:lightgray, showsegments=false, floor_offset=-20, ceiling_offset=10)
        plotmesh!(ax3, mm;color=Za, showsegments=false, floor_offset=-20, ceiling_offset=10)
    end
    Colorbar(lg[1,7], colorrange=(extrema(filter(isfinite, Za))), label="Firing rate [Hz]")
    if embeddim(mm) == 2
    for ax in [ax1, ax2, ax3]
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.xticksvisible = false
        ax.yticksvisible = false
        ax.bottomspinevisible = false
        ax.leftspinevisible = false
        end
        if !isempty(ylabel)
            ax1.ylabel = ylabel
        end
    end
    colsize!(lg, 1, Relative(0.1))
    ax1,ax2,ax3
end

"""
Plot a summary of the spatial selectivity across all cells

One panel showing the distribution of SIC scores
Some example cells. One near the top, i.e. the most spatially selective, one near the median and one near the bottom?
"""
function plot_spatial_summary!(lg,celldirs::Union{Vector{String},Nothing}=nothing;redo=false, kwargs...)
    if celldirs === nothing
        celldirs = open("/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/cell_list.txt") do fid
                    readlines(fid)
       end
    end
    h = UInt32(0) 
    for c in celldirs
        h = CRC32c.crc32c(c,h)
    end
    hs = string(h,base=16)
    fname = "spatial_summary_data_$(hs).jld2"
    if !redo && isfile(fname)
        res,sic,mean_fr = JLD2.load(fname, "res","sic","mean_fr")
    else
        res = issignificant(SpatialInformationContent, celldirs;skip_error=true, kwargs...)
        sic = get_sic(SpatialInformationContent, celldirs;skip_error=true, kwargs...)
        mean_fr = process_dirs((;kwargs...)->get_mean_firing_rate(JointMap(;kwargs...)), celldirs;skip_error=true)
        # reformat to array
        mean_fr = [get(mean_fr, c, NaN) for c in celldirs]
        JLD2.save(fname, Dict("res"=>res, "sic"=>sic,"mean_fr"=>mean_fr))
    end

    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    ax1 = Axis(lg1[1,1])
    _colors = parse.(Colorant, [:grey, :red])
    cc = fill(_colors[1], length(sic))
    cc[res] .= _colors[2] 
    xx = rand(length(cc))
    yy = sic
    scatter!(ax1, mean_fr, yy, color=cc)
    Legend(lg1[1,2], [MarkerElement(marker=:circle, color=q) for q in _colors],["Non-sig","sig"],tellwidth=false, tellheight=false,
           valign=:top, halign=:right, framevisible=true, padding=(10.0, 10.0, 10.0, 10.0))
    #rainclouds!(ax1, xx, yy)
    #ax1.xticks = (1:2, ["Non-sig", "sig"])
    ax1.ylabel = "SIC"
    ax1.xlabel = "Mean firing rate [Hz]"
    ax1.xticklabelsvisible = true 
    ax1.xticksvisible = true 
    ax1.bottomspinevisible = true 
    ax2 = Axis(lg1[1,2])
    linkyaxes!(ax1,ax2)
    density!(ax2, yy[res.==false], direction=:y, color=_colors[1])
    density!(ax2, yy[res.==true], direction=:y, color=_colors[2],alpha=0.8)
    ax2.yticklabelsvisible = false
    ax2.xticklabelsvisible = false
    ax2.xticksvisible = false
    ax3 = Axis(lg1[1,3],aspect=1)
    pie!(ax3, [sum(res.==true), sum(res.==false)], color=[:gray, :red])
    aa = round(100*sum(res)/length(res), sigdigits=3)
    bb = round(100-aa, sigdigits=3)
    text!(ax3, [0.5, 0.6],[0.75, 0.25], text=["$aa%","$bb%"], color=:white,space=:relative)
    hidedecorations!(ax3)
    ax3.bottomspinevisible = false
    ax3.leftspinevisible = false

    place_idx = findall(res.==true)
    non_place_idx = findall(res.==false)
    # TODO: Include an example
    lg2 = GridLayout(lg[2,1])
    Label(lg2[1,1,TopLeft()], "B")
    cidx1 = place_idx[argmax(sic[place_idx])]
    plot_place_map!(lg2, celldirs[cidx1])
    Label(lg2[1,0], "Cell 1", rotation=π/2, tellheight=false,color=:blue)

    lg3 = GridLayout(lg[3,1])
    cidx2 = non_place_idx[argmax(sic[non_place_idx])]
    plot_place_map!(lg3, celldirs[cidx2])
    Label(lg3[1,0], "Cell 2", rotation=π/2, tellheight=false,color=:green)

    lg4 = GridLayout(lg[4,1])
    # median
    cidx3 = place_idx[sortperm(sic[place_idx])[div(length(place_idx),2)]]
    plot_place_map!(lg4, celldirs[cidx3];ylabel="Cell 3")
    Label(lg4[1,0], "Cell 3", rotation=π/2, tellheight=false,color=:purple)

    # TODO: Indicate these two cells in the scatter plot
    cidx = [cidx1, cidx2, cidx3]
    scatter!(ax1, mean_fr[cidx], yy[cidx], markersize=10px, strokecolor=[:blue,:green,:purple], color=cc[cidx], strokewidth=2.0)
    annotation!(ax1, mean_fr[cidx] .+ [15, 15, 0.0], yy[cidx] .+ [10.0, 0.0, -15], mean_fr[cidx], yy[cidx]; text=["1","2","3"],textcolor=[:blue, :green, :purple]) 
end

function plot_spatial_summary(celldirs::Union{Vector{String},Nothing}=nothing;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(1035, 1034))
        lg = GridLayout(fig[1,1])
        plot_spatial_summary!(lg, celldirs;kwargs...)
        fig
    end
end

function plot_response_fields(::Type{T}, args...;kwargs...) where T <: AbstractResponseFields
    ncells = length(args[1])
    if length(args) == 2
        width = 300*length(args[2])
    else
        width = 600
    end
    height = 400*ncells
    with_theme(plot_theme) do
        fig = Figure(size=(width,height))
        lg = GridLayout(fig[1,1])
        plot_response_fields!(lg, T, args...;kwargs...)
        fig
    end
end

function plot_response_fields!(lg, ::Type{T}, celldirs::Vector{String},args...;kwargs...) where T <: AbstractResponseFields
    lg = [GridLayout(lg[i,1]) for i in 1:length(celldirs)]
    labels = range('A', length=length(celldirs), step=1)
    for (_lg, celldir,ll) in zip(lg,celldirs,labels)
        plot_response_fields!(_lg, T, celldir,args...;kwargs...)
        Label(_lg[1,1,TopLeft()], string(ll))
    end
end

function plot_response_fields!(lg, ::Type{T}, celldir::String;nshuffles=10_000, nrefinements=(p=3,g=2), kwargs...) where T <: AbstractResponseFields
    if T <: SpatialResponseFields
        TM = SpatialInformationContent
    else
        TM = GazeInformationContentent
    end
    mm = get_mesh(T,nrefinements)
    rf,sic = cd(celldir) do 
        rf = get_response_fields(T, nshuffles;nrefinements=nrefinements,kwargs...)
        sic = compute_skaggs_sic(TM,nshuffles;nrefinements=nrefinements,kwargs...)
        rf,sic
    end
    Z = zeros(nelements(mm))
    Z[rf.binidx] .= 1.0
    
    if embeddim(mm) == 3
        lscene = LScene(lg[1,1])
        plotmesh!(lscene, mm;showsegments=true, color=Z)
    else
        ax = Axis(lg[1,1])
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        viz!(ax, mm;showsegments=true, color=Z)
    end
    # show distribution of sic
    ax2 = Axis(lg[2,1])
    boxplot!(ax2, fill(1.0, length(sic.sic)), sic.sic;show_outliers=false, show_notch=true,color=:gray,width=0.8)
    scatter!(ax2, [1.0], [sic.sic0], color=:red)
    ax2.xticksvisible = false
    ax2.xticklabelsvisible = false
    ax2.bottomspinevisible = false
    ylabelvisible = get(kwargs, :ylabelvisible, false)
    if ylabelvisible
        ax2.ylabel = "SIC"
    end
end

function plot_response_fields!(lg, ::Type{T}, celldir::String, nrefinements::Vector{T2};kwargs...) where T <: AbstractResponseFields where T2 <: @NamedTuple{p::Int64, g::Int64}
    lgs = [GridLayout(lg[1,i]) for i in 1:length(nrefinements)]
    for (ii,(_lg, _nrefinements)) in enumerate(zip(lgs, nrefinements))
        ylabelvisible = ii == 1 
        plot_response_fields!(_lg, T, celldir;nrefinements=_nrefinements,ylabelvisible=ylabelvisible, kwargs...)
    end
end

function plot_map(axes, ii, mp::Type{T}, celldirs::Vector{String};nrefinements=(p=3,g=2), kwargs...) where T <: AbstractMap
    # axes[1] for map axes[2] for colorbar
    mm = get_mesh(T;nrefinements=nrefinements)
    Z = lift(ii) do i
        jm = cd(celldirs[i]) do
            JointMap(;nrefinements=nrefinements, kwargs...)
        end
        spm = T(jm,mm)
        Z = spm.weight./spm.occupancy
        zmin,zmax = extrema(filter(isfinite, Z))
        axes[2].colorrange[] = (zmin,zmax)
        Z
    end
    if embeddim(mm) == 2
        viz!(axes[1],mm;showsegments=true,color=Z)
    else
        plotmesh!(axes[1],mm;showsegments=true, color=Z, floor_offset=-20, ceiling_offset=10)
    end
    length(celldirs)
end

plot_spatial_map(axes, ii, celldirs;kwargs...) = plot_map(axes, ii, SpatialMapNew, celldirs;kwargs...)
plot_gaze_map(axes, ii, celldirs;kwargs...) = plot_map(axes, ii, ViewMapNew, celldirs;kwargs...)

function plot_spatial_map(lg::GridLayout)
    ax = Axis(lg[1,1])
    axc = Colorbar(lg[1,2], colorrange=(0.0, 0.0), label="Firing rate [Hz]")
    ax,axc
end


function plot_gaze_map(lg::GridLayout)
    ax = LScene(lg[1,1])
    axc = Colorbar(lg[1,2], colorrange=(0.0, 0.0), label="Firing rate [Hz]")
    ax,axc
end

function plot_conjunctions(λ_covered, λ_not_covered, mm)
    with_theme(plot_theme) do
        fig = Figure()
        lscene1 = LScene(fig[2,1])
        lscene2 = LScene(fig[2,2])

        x1min,x1max = extrema(filter(isfinite,λ_covered))
        x2min,x2max = extrema(filter(isfinite,λ_not_covered))
        xmin = min(x1min, x2min)
        xmax = max(x1max, x2max)
        plotmesh!(lscene1, mm;color=:lightgray, floor_offset=-20, ceiling_offset=10)
        plotmesh!(lscene1, mm;color=λ_covered, floor_offset=-20, ceiling_offset=10,colorrange=(xmin,xmax))

        plotmesh!(lscene2, mm;color=:lightgray, floor_offset=-20, ceiling_offset=10)
        plotmesh!(lscene2, mm;color=λ_not_covered, floor_offset=-20, ceiling_offset=10, colorrange=(xmin,xmax))
        Colorbar(fig[2,3], colorrange=(xmin,xmax), label="Firing Rate [Hz]")
        Label(fig[1,1], "In field", tellwidth=false)
        Label(fig[1,2], "Out-of field", tellwidth=false)
        fig
    end
end

function plot_fields(::Type{T}, celldir::String;colormap=:binary, kwargs...) where T <: AbstractResponseFields
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    jm,rf = cd(celldir) do
        jm = JointMap(;kwargs...)
        rf = get_response_fields(T, 10_000;nrefinements=nrefinements, smooth=true, smoothing_method=:laplace,
                                                                α=0.1, niter=100)
        jm,rf
    end
    mm = get_mesh(T,nrefinements)
    # cluster fields based on similarity
    clusters = merge_fields(mm, rf.binidx)
    if length(clusters) <= 7
        wcolors = Makie.wong_colors()
    else
        wcolors = to_colormap(:Paired_12)
    end
    cluster_colors = fill(wcolors[1], length(rf.binidx))
    for (ii,cluster) in enumerate(clusters)
        cluster_colors[cluster] .= wcolors[ii] 
    end
    TM = maptype(T)
    spm = TM(jm, mm)
    sml = SmoothedMap(spm;method=:laplace, α=0.1, niter=100)

    with_theme(plot_theme) do
        fig = Figure(size=(800,400))
        if embeddim(mm) == 2
            ax1 = Axis(fig[1,1])
            ax2 = Axis(fig[1,3])
              for ax in [ax1,ax2]
                hidedecorations!(ax)
                ax.topspinevisible = true
                ax.rightspinevisible = true
            end
        else
            ax1 = LScene(fig[1,1],show_axis=false)
            ax2 = LScene(fig[1,3], show_axis=false)
        end
        if embeddim(mm) == 3
            plotmesh!(ax1, mm;color=get_rate_map(spm),colormap=colormap,floor_offset=-20, ceiling_offset=10,showsegments=true,segmentcolor=:lightgray)
            plotmesh!(ax2, mm;color=get_rate_map(sml),colormap=colormap, floor_offset=-20, ceiling_offset=10,showsegments=true,segmentcolor=:lightgray)
            m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
            ceiling_idx = findall(in(m_ceiling.inds).(rf.binidx))
            floor_idx = findall(in(m_floor.inds).(rf.binidx))
            centroids = centroid.(mm[rf.binidx])
            centroids[ceiling_idx] .= Translate(0.0, 0.0, 10).(centroids[ceiling_idx])
            centroids[floor_idx] .= Translate(0.0, 0.0, -20).(centroids[floor_idx])
            viz!(ax2, centroids;color=cluster_colors,pointsize=7.5)
        else
            viz!(ax1, mm;color=get_rate_map(spm),colormap=colormap)
            viz!(ax2, mm;color=get_rate_map(sml),colormap=colormap)
            viz!(ax2, centroid.(mm[rf.binidx]),color=cluster_colors)
        end
        Colorbar(fig[1,2], colorrange=extrema(filter(isfinite, get_rate_map(spm))), label="Firing rate [Hz]",ticksvisible=true,colormap=colormap)
        Colorbar(fig[1,4], colorrange=extrema(filter(isfinite, get_rate_map(sml))), label="Firing rate [Hz]", ticksvisible=true, colormap=colormap)
        link_cameras_lscene(fig)
        fig
    end
end

function plot_laplace_adaptive_comparison(celldir::String;laplace_params=(α=0.1,niter=100), adaptive_params=(;α=1000.0^2))
    spm = cd(celldir) do
        SpatialMapNew()
    end
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));
    sml = SmoothedMap(spm;method=:laplace, laplace_params...)
    sma = SmoothedMap(spm;method=:adaptive, adaptive_params...)
    with_theme(plot_theme) do
        fig = Figure(size=(800,300))

        axes = [Axis(fig[1,2*(i-1)+1],aspect=1, rightspinevisible=true, topspinevisible=true) for i in 1:3]
        hidedecorations!.(axes)
        axes[1].title = L"Raw"
        Z = get_rate_map(spm)
        viz!(axes[1],m_floor, color=Z)
        Colorbar(fig[1,2], colorrange=extrema(filter(isfinite, Z)))
        axes[2].title = L"Adaptive, \alpha=%$(adaptive_params[:α])"
        Za = get_rate_map(sma)
        viz!(axes[2], m_floor, color=Za)
        Colorbar(fig[1,4], colorrange=extrema(filter(isfinite, Za)))
        axes[3].title =L"Laplace, \alpha=%$(laplace_params[:α]), niter=%$(laplace_params[:niter])"
        Zl = get_rate_map(sml)
        viz!(axes[3],m_floor, color=Zl)
        Colorbar(fig[1,6], colorrange=extrema(filter(isfinite, Zl)))
        fig
    end
end

function plot_maps(jm::JointMap;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_maps!(lg, jm;kwargs...)
        fig
    end
end

function plot_maps!(lg, jm::JointMap;smooth=true, smoothing_method=:laplace, α=0.1, niter=100)

    mm = get_maze_mesh(;nrefinements=jm.dims[1])
    m_floor = floor_topology3(;nrefinements=jm.dims[2])
    # view map
    vm = ViewMapNew(jm, mm)
    # place map
    spm = SpatialMapNew(jm, Shadow("xy")(m_floor))
    if smooth
        vml = SmoothedMap(vm;method=smoothing_method, α=α, niter=niter)
        sml= SmoothedMap(spm;method=smoothing_method, α=α, niter=niter)
    else
        vml = vm
        sml = spm 
    end
    lscene = LScene(lg[1,1], show_axis=false)
    λv = get_rate_map(vml)
    λsp = get_rate_map(sml)
    plotmesh!(lscene,mm;color=:lightgray, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene,mm;color=λv, ceiling_offset=10, floor_offset=-15)
    # offset the floor
    m_floor = Translate(0.0, 0.0, -30)(m_floor)
    viz!(lscene, m_floor;color=:lightgray)
    viz!(lscene, m_floor;color=λsp)
end