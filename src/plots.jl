using CairoMakie


function plot_map(::Type{T}, celldir::String;kwargs...) where T <: AbstractInformationContent
    with_theme(merge(plot_theme, theme_latexfonts())) do
        fig = Figure(size=(1000,500))
        lg = GridLayout(fig[1,1])
        plot_map!(lg, T, celldir;kwargs...)
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
    offset = 0
    axes = Any[]
    for (ii, (_sm, _sic,sparams)) in enumerate(zip([sm, smg, sma],[sic, sicg,sica],["raw";smooth_params]))
        Z = get_rate_map(_sm)
        if embeddim(mm) == 2
            ax = Axis(lg[1,offset+2*(ii-1)+1])
            ax.xticklabelsvisible = false
            ax.yticklabelsvisible = false
            ax.xticksvisible = false
            ax.yticksvisible = false
            ax.bottomspinevisible = false
            ax.leftspinevisible = false
            if !isempty(ylabel)
                ax1.ylabel = ylabel
            end
            viz!(ax, mm;showsegments=false, color=:lightgray)
            viz!(ax, mm;showsegments=false, color=Z)
            ax.title = strip(string(sparams),['(',')']) 
        else
            ax = LScene(lg[1,offset+2*(ii-1)+1],show_axis=false)
            plotmesh!(ax, mm;color=:lightgray, showsegments=false, floor_offset=-20, ceiling_offset=10)
            plotmesh!(ax, mm;color=Z, showsegments=false, floor_offset=-20, ceiling_offset=10)
        end
        push!(axes, ax)
        Colorbar(lg[1,offset+2*(ii-1)+2], colorrange=(extrema(filter(isfinite, Z))),label="Firing rate [Hz]")
        # plot sic distribution below
        ax2 = Axis(lg[2, offset+2*(ii-1)+1],ytickformat=value->["$(round(v*100, sigdigits=2))" for v in value])
        Label(lg[2, offset+2*(ii-1)+1, TopLeft()], L"\times 10^{-2}")
        ax2.xticklabelsvisible = false
        ax2.xticksvisible = false
        ax2.bottomspinevisible = false
        boxplot!(ax2,  fill(1.0, length(_sic.sic)), _sic.sic;show_outliers=false, show_notch=true,color=:gray)
        hlines!(ax2, percentile(_sic.sic,95), linestyle=:dot)
        scatter!(ax2, [1.0],[_sic.sic0],color=:red)
        if ii == 1
            ax2.ylabel = "SIC"
        end
    end
    rowsize!(lg, 1, Relative(0.6))
    axes
end

"""
Plot a summary of the spatial selectivity across all cells

One panel showing the distribution of SIC scores
Some example cells. One near the top, i.e. the most spatially selective, one near the median and one near the bottom?
"""
function plot_summary!(lg,::Type{T}, celldirs::Union{Vector{String},Nothing}=nothing;redo=false, kwargs...) where T <: AbstractInformationContent
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
    if T <: SpatialInformationContent
        fname = "spatial_summary_data_$(hs).jld2"
    else
        fname = "gaze_summary_data_$(hs).jld2"
    end
    if !redo && isfile(fname)
        res,sic,mean_fr = JLD2.load(fname, "res","sic","mean_fr")
    else
        res = issignificant(T, celldirs;skip_error=true, kwargs...)
        sic = get_sic(T, celldirs;skip_error=true, kwargs...)
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
    pie!(ax3, [sum(res.==true), sum(res.==false)], color=[:red, :gray])
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
    plot_map!(lg2, T, celldirs[cidx1])
    Label(lg2[1,0], "Cell 1", rotation=π/2, tellheight=false,color=:blue)

    lg3 = GridLayout(lg[3,1])
    cidx2 = non_place_idx[argmax(sic[non_place_idx])]
    plot_map!(lg3, T, celldirs[cidx2])
    Label(lg3[1,0], "Cell 2", rotation=π/2, tellheight=false,color=:green)

    lg4 = GridLayout(lg[4,1])
    # median
    cidx3 = place_idx[sortperm(sic[place_idx])[div(length(place_idx),2)]]
    plot_map!(lg4, T, celldirs[cidx3];ylabel="Cell 3")
    Label(lg4[1,0], "Cell 3", rotation=π/2, tellheight=false,color=:purple)

    # TODO: Indicate these two cells in the scatter plot
    cidx = [cidx1, cidx2, cidx3]
    @show celldirs[cidx]
    scatter!(ax1, mean_fr[cidx], yy[cidx], markersize=10px, strokecolor=[:blue,:green,:purple], color=cc[cidx], strokewidth=2.0)
    annotation!(ax1, mean_fr[cidx] .+ [15, 15, 0.0], yy[cidx] .+ [10.0, 0.0, -15], mean_fr[cidx], yy[cidx]; text=["1","2","3"],textcolor=[:blue, :green, :purple]) 
end

function plot_summary(::Type{T}, celldirs::Union{Vector{String},Nothing}=nothing;kwargs...) where T <: AbstractInformationContent
    with_theme(plot_theme) do
        fig = Figure(size=(1035, 1034))
        lg = GridLayout(fig[1,1])
        plot_summary!(lg, T, celldirs;kwargs...)
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
        TM = GazeInformationContent
    end
    mm = get_mesh(T,nrefinements)
    rf,sic = cd(celldir) do 
        rf = get_response_fields(T, nshuffles;nrefinements=nrefinements,kwargs...)
        sic = compute_skaggs_sic(TM,nshuffles;nrefinements=nrefinements,kwargs...)
        rf,sic
    end
    Z = zeros(nelements(mm))
    Z[rf.binidx] .= 1.0
    @show rf.binidx 
    if embeddim(mm) == 3
        lscene = LScene(lg[1,1])
        plotmesh!(lscene, mm;showsegments=true, color=Z,floor_offset=-20, ceiling_offset=10)
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
    hlines!(percentile(sic.sic, 95), linestyle=:dot)
    ax2.xticksvisible = false
    ax2.xticklabelsvisible = false
    ax2.bottomspinevisible = false
    ylabelvisible = get(kwargs, :ylabelvisible, false)
    if ylabelvisible
        ax2.ylabel = "SIC"
    end
    rowsize!(lg, 1, Relative(0.7))
    # return the elements
    mm[rf.binidx], issignificant(sic)
end

function plot_response_fields!(lg, ::Type{T}, celldir::String, nrefinements::Vector{T2};kwargs...) where T <: AbstractResponseFields where T2 <: @NamedTuple{p::Int64, g::Int64}
    lgs = [GridLayout(lg[1,i]) for i in 1:length(nrefinements)]
    qq = Any[]
    for (ii,(_lg, _nrefinements)) in enumerate(zip(lgs, nrefinements))
        ylabelvisible = ii == 1 
        _qq,_sig = plot_response_fields!(_lg, T, celldir;nrefinements=_nrefinements,ylabelvisible=ylabelvisible, kwargs...)
        if _sig
            append!(qq, _qq)
        end
    end
    if !isempty(qq)
        mm = get_mesh(T, (p=0,g=0))
        qq = convert(Vector{typeof(qq[1])}, qq)
        # consolidate the fields
        qqf = merge_fields(qq)
        if embeddim(mm) == 2
            ax = Axis(lg[1, length(lgs)+1], aspect=1)
            ax.xticklabelsvisible = false
            ax.yticklabelsvisible = false
        else
            ax = LScene(lg[1, length(lgs)+1])
        end
        colsize!(lg, length(lgs)+1, Relative(0.3))
        # figure out the colors
        if 3 <= length(qqf) <= 12
            colors = to_colormap(Symbol("Paired_$(length(qqf))"))
        else
            colors = resample_cmap(:Paired_12, length(qqf))
        end
        if embeddim(mm) == 2
            bb = Meshes.boundingbox(mm)
            viz!(ax, bb, color=:black)
            viz!(ax, mm;showsegments=false, color=:white)
            viz!(ax, boundary(bb), color=:black)
            viz!(ax, qqf, color=colors,linewidth=3.0)
            ax.title = "# fields: $(length(qqf))"
        else
            # TODO: Offset patches on the floor and ceiling
            m_floor, m_ceiling, m_middle = get_floor_and_ceiling(mm)
            floor_idx = findall(Meshes.intersects(m_floor), qqf)
            ceil_idx = findall(Meshes.intersects(m_ceiling), qqf)
            plotmesh!(ax, mm;showsegments=false, alpha=1.0, floor_offset=0, ceiling_offset=0,hide_ceiling=true,hide_floor=true, alphasegments=0.5)
            #qqf[floor_idx] = Translate(0.0, 0.0, -20).(qqf[floor_idx])
            #qqf[ceil_idx] = Translate(0.0, 0.0, 10).(qqf[ceil_idx])
            viz!(ax, qqf, color=colors, linewidth=3.0)
        end
    end
end

function plot_response_fields(ax, ii, rfs::Vector{T},maps::Union{Nothing, Vector{<:AbstractMap}}=nothing;showmaze=true, kwargs...) where T <: AbstractResponseFields
    nrefinements = first(rfs).args[:nrefinements]
    mm = get_mesh(T,nrefinements)
    nd = embeddim(mm)
    cm = to_colormap(:Paired_12)
    points = centroid.(mm)
    colors = Observable(fill(parse(Colorant, :lightgray), length(points)))
    alpha = Observable(fill(0.0, length(points)))
    Z = Observable(zeros(nelements(mm)))
    on(ii)  do _ii
        _rf = rfs[_ii]
        clusters = merge_fields(_rf)
        _colors = colors[]
        _alpha = alpha[]
        fill!(_colors, parse(Colorant, :lightgray))
        fill!(_alpha, 0.0)
        if !isempty(clusters)
             for (ii,cluster) in enumerate(clusters)
                _colors[_rf.binidx[cluster]] .= cm[ii]
                _alpha[_rf.binidx[cluster]] .= 1.0
             end
        end
        colors[] .= _colors
        alpha[] .= _alpha
        notify(colors)
        notify(alpha)
        if maps !== nothing
            Z[] = get_rate_map(maps[_ii])
        end
    end
    if maps !== nothing
        if embeddim(mm) == 2
            viz!(ax, mm;color=Z,colormap=:binary)
        else
            plotmesh!(ax, mm;color=Z)
        end
    else
        if showmaze
            viz!(ax, mm;color=:lightgray)
        end
    end
    viz!(ax, points;color=colors, pointsize=10.0,alpha=alpha)
    return length(rfs)
end

function plot_response_fields(lg::GridLayout)
    ax = Axis(lg[1,1])
    ax
end

function plot_map(axes, ii, mp::Type{T}, celldirs::Vector{String};nrefinements=(p=3,g=2), kwargs...) where T <: AbstractMap
    # axes[1] for map axes[2] for colorbar
    mm = get_mesh(T;nrefinements=nrefinements)
    Z = lift(ii) do i
        jm = cd(celldirs[i]) do
            JointMap(;nrefinements=nrefinements, kwargs...)
        end
        spm = T(jm,mm)
        if get(kwargs, :smooth, false)
            spm = SmoothedMap(spm;kwargs...)
            Z = spm.weight./spm.occupancy
            Z[spm.unvisited] .= NaN
        else
            Z = spm.weight./spm.occupancy
        end
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

function plot_maps!(lg, jm::JointMap;sic_spatial::Union{SpatialInformationContent,Nothing}=nothing,
                                      sic_gaze::Union{GazeInformationContent,Nothing}=nothing,
                                      rf_spatial::Union{SpatialResponseFields,Nothing}=nothing,
                                      rf_gaze::Union{GazeResponseFields, Nothing}=nothing,
                                      smooth=true, smoothing_method=:laplace, α=0.1, niter=100,colormap=:binary)

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
    plotmesh!(lscene,mm;color=:lightblue, ceiling_offset=10, floor_offset=-15)
    plotmesh!(lscene,mm;color=λv, ceiling_offset=10, floor_offset=-15,colormap=colormap)
    if rf_gaze !== nothing
        plot_response_fields!(lscene, rf_gaze;floor_offset=-15, ceiling_offset=10)
    end
    # offset the floor
    m_floor = Translate(0.0, 0.0, -30)(m_floor)
    viz!(lscene, m_floor;color=:lightblue)
    viz!(lscene, m_floor;color=λsp,colormap=colormap)

    if rf_spatial !== nothing
        plot_response_fields!(lscene, rf_spatial;offset=-30)
    end
    # show colorbar
    lg2 = GridLayout(lg[1,2])
    Colorbar(lg2[1,1], colorrange=extrema(filter(isfinite, λv)),label="Firing rate [Hz]",colormap=colormap)
    Colorbar(lg2[2,1], colorrange=extrema(filter(isfinite, λsp)), label="Firing rate [Hz]",colormap=colormap)
    rowsize!(lg2, 1, Relative(0.7))
    if sic_spatial !== nothing || sic_gaze !== nothing
        lg3 = GridLayout(lg[1,3])
        colsize!(lg, 3, 75)
        ax1 = Axis(lg3[1,1])
        boxplot!(ax1, fill(1.0, length(sic_spatial.sic)), sic_spatial.sic, show_outliers=false, show_notch=true, color=:darkgray)
        hlines!(ax1, sic_spatial.sic0, linestyle=:dot, color=:black)

        ax2 = Axis(lg3[2,1])
        boxplot!(ax2, fill(1.0, length(sic_gaze.sic)), sic_gaze.sic, show_outliers=false, show_notch=true, color=:darkgray)
        hlines!(ax2, sic_gaze.sic0, linestyle=:dot, color=:black)
        for (ll,_ax) in zip(["View","Place"], [ax1, ax2])
            _ax.bottomspinevisible = false
            _ax.xticksvisible = false
            _ax.xticklabelsvisible = false
            _ax.yaxisposition = :right
            _ax.leftspinevisible = false
            _ax.rightspinevisible = true
            _ax.ylabel = "SIC $ll"
        end
    end

end


function plot_trajectories_with_head_direction(pos::Vector{<:Matrix{<:Real}}, head_direction::Vector{<:Vector{<:Real}};kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_trajectories_with_head_direction!(lg, pos, head_direction;kwargs...)
        fig
    end
end

function plot_trajectories_with_head_direction!(lg, pos::Vector{<:Matrix{<:Real}}, head_direction::Vector{<:Vector{<:Real}};placefieldidx::Union{Nothing, Vector{Int64}}=nothing)
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));
    ax = Axis(lg[1,1], aspect=1) 
    hidedecorations!(ax)
    ax.bottomspinevisible = false
    ax.leftspinevisible = false
    viz!(ax, m_floor;color=:lightgray)
    for (_pos, hd) in zip(pos, head_direction)
        lines!(ax, Point2f.(eachcol(_pos[1:2,:])), color=hd, colormap=:phase, linewidth=2.0)
    end
    if placefieldidx !== nothing
        bb = find_boundary(m_floor[placefieldidx])
        viz!(ax, bb, color=:black)
    end
    axp = Axis(lg[1,1], width=Relative(0.2), height=Relative(0.2), halign=0.85, valign=0.05, aspect=1)
    hidedecorations!(axp)
    axp.backgroundcolor = RGB(0.8, 0.8, 0.8) 
    axp.topspinevisible = true 
    axp.rightspinevisible = true 
    axp.yaxisposition = :right
    axp.ylabel = "Head direction"
    axp.ylabelvisible = true
    θ = range(0.0, stop=2π, length=100)
    lines!(axp, cos.(θ), sin.(θ), color=θ, colormap=:phase, linewidth=4.0)
    # make sure the inset is on top
    translate!(axp.blockscene, 0, 0, 150)
end

function plot_spatial_responses_comparison(celldir)
    kwargs1 = (refinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0,trial_start=2, pv_threshold=0.001)
    kwargs2 = (refinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)

    with_theme(plot_theme)  do
        fig = Figure(size=(900,1200))
        lg1 = GridLayout(fig[1,1])
        Label(lg1[1,1,TopLeft()], "A")
        plot_spatial_responses!(lg1, celldir;kwargs1...)
        lg2 = GridLayout(fig[2,1])
        Label(lg2[1,1,TopLeft()], "B")
        plot_spatial_responses!(lg2, celldir;kwargs2...)
        fig
    end
end
function plot_spatial_responses(celldir::String;kwargs...)
     with_theme(plot_theme) do
        fig = Figure(size=(1000,700))
        lg = GridLayout(fig[1,1]) 
        plot_spatial_responses!(lg, celldir;kwargs...)
        fig
    end
end

function plot_spatial_responses!(lg, celldir::String;kwargs...)
    
        jm,rf = cd(celldir) do 
            jm = JointMap(;kwargs...)
            rf = get_response_fields(SpatialResponseFields, 10_000;kwargs...)
            jm, rf
        end
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=3))
        spm = SpatialMapNew(jm, m_floor)
        colormap = get(kwargs, :colormap, :rain)
        # smoothing
        sml = SmoothedMap(spm;method=:laplace, α=0.1, niter=100)
        axes = [Axis(lg[1,i], aspect=1, xticklabelsvisible=false, yticklabelsvisible=false,
                                xticksvisible=false, yticksvisible=false) for i in 1:3]
        viz!(axes[1], m_floor, color=:lightgray,showsegments=false)
        Z = spm.occupancy
        viz!(axes[1], m_floor, color=Z,colormap=colormap)
        Colorbar(lg[2,1], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Occupancy [s]", colormap=colormap)

        Z = spm.weight
        viz!(axes[2], m_floor, color=:lightgray,showsegments=false)
        viz!(axes[2], m_floor, color=Z, colormap=colormap)
        Colorbar(lg[2,2], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Spike count", colormap=colormap)

        Z = get_rate_map(spm)
        viz!(axes[3], m_floor, color=:lightgray,showsegments=false)
        viz!(axes[3], m_floor, color=Z, colormap=colormap)
        Colorbar(lg[2,3], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Firing rate [Hz]", colormap=colormap)


        # plot the equivalent for the smoothed responses
        axes2 = [Axis(lg[3,i], aspect=1, xticklabelsvisible=false, yticklabelsvisible=false,
                                xticksvisible=false, yticksvisible=false) for i in 1:3] 

        Z = sml.occupancy
        viz!(axes2[1], m_floor, color=:lightgray,showsegments=false)
        viz!(axes2[1], m_floor, color=Z,colormap=colormap)
        Colorbar(lg[4,1], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Occupancy [s]", colormap=colormap)

        Z = sml.weight
        viz!(axes2[2], m_floor, color=:lightgray,showsegments=false)
        viz!(axes2[2], m_floor, color=Z, colormap=colormap)
        Colorbar(lg[4,2], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Spike count", colormap=colormap)

        Z = get_rate_map(sml)

        #viz!(axes[4], m_floor, color=:lightgray,showsegments=false)
        #viz!(axes[4], m_floor, color=Z, colormap=colormap)
        #Colorbar(fig[2,4], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Firing rate [Hz]", colormap=colormap)
        lg1 = GridLayout(lg[3:4,3])        
        plot_response_fields!(lg1, rf;filter_spurious=true, colorbar_below=true, colormap=:rain)
end

function plot_view_response_comparison(celldir::String)
    kwargs1 = (refinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0,trial_start=2, pv_threshold=0.001)
    kwargs2 = (refinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)

    with_theme(plot_theme)  do
        fig = Figure(size=(900,1200))
        lg1 = GridLayout(fig[1,1])
        Label(lg1[1,1,TopLeft()], "A")
        plot_view_responses!(lg1, celldir;kwargs1...)
        lg2 = GridLayout(fig[2,1])
        Label(lg2[1,1,TopLeft()], "B")
        plot_view_responses!(lg2, celldir;kwargs2...)
        fig
    end
end

function plot_view_responses(celldir::String;kwargs...)
     with_theme(plot_theme) do
        fig = Figure(size=(1000,800))
        lg = GridLayout(fig[1,1]) 
        plot_view_responses!(lg, celldir;kwargs...)
        fig
    end
end

function plot_view_responses!(lg, celldir::String;do_animate=false, fname::Union{String, Nothing}=nothing, kwargs...)
    jm,rf = cd(celldir) do 
        jm = JointMap(;kwargs...)
        rf = get_response_fields(GazeResponseFields, 10_000;kwargs...)
        jm, rf
    end
    mm = get_maze_mesh(;nrefinements=2)
    spm = ViewMapNew(jm, mm)
    colormap = get(kwargs, :colormap, :rain)
    floor_offset = get(kwargs, :floor_offset,-15)
    ceiling_offset = get(kwargs, :ceiling_offset,15)
    # smoothing
    sml = SmoothedMap(spm;method=:laplace, α=0.1, niter=100)
    lscene1 = LScene(lg[1,1], show_axis=false)
    plotmesh!(lscene1, mm, color=:lightgray,showsegments=true,floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Z = spm.occupancy
    plotmesh!(lscene1, mm, color=Z,colormap=colormap, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Colorbar(lg[1,2], colorrange=extrema(filter(isfinite, Z)), label="Occupancy [s]", colormap=colormap)

    lscene2 = LScene(lg[1,3], show_axis=false)
    Z = spm.weight
    plotmesh!(lscene2, mm, color=:lightgray,showsegments=true, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    plotmesh!(lscene2, mm, color=Z, colormap=colormap, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Colorbar(lg[1,4], colorrange=extrema(filter(isfinite, Z)), label="Spike count", colormap=colormap)

    lscene3 = LScene(lg[1,5],show_axis=false)
    Z = get_rate_map(spm)
    plotmesh!(lscene3, mm, color=:lightgray,showsegments=true,floor_offset=floor_offset, ceiling_offset=ceiling_offset )
    plotmesh!(lscene3, mm, color=Z, colormap=colormap,floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Colorbar(lg[1,6], colorrange=extrema(filter(isfinite, Z)), label="Firing rate [Hz]", colormap=colormap)


    # plot the equivalent for the smoothed responses
    lscene4 = LScene(lg[2,1], show_axis=false) 
    Z = sml.occupancy
    plotmesh!(lscene4, mm, color=:lightgray,showsegments=true,floor_offset=floor_offset, ceiling_offset=ceiling_offset )
    plotmesh!(lscene4, mm, color=Z,colormap=colormap, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Colorbar(lg[2,2], colorrange=extrema(filter(isfinite, Z)), label="Occupancy [s]", colormap=colormap)

    lscene5 = LScene(lg[2,3],show_axis=false)
    Z = sml.weight
    plotmesh!(lscene5, mm, color=:lightgray,showsegments=true, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    plotmesh!(lscene5, mm, color=Z, colormap=colormap, floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    Colorbar(lg[2,4], colorrange=extrema(filter(isfinite, Z)), label="Spike count", colormap=colormap)

    Z = get_rate_map(sml)

    #viz!(axes[4], m_floor, color=:lightgray,showsegments=false)
    #viz!(axes[4], m_floor, color=Z, colormap=colormap)
    #Colorbar(fig[2,4], colorrange=extrema(filter(isfinite, Z)), vertical=false, flipaxis=false, label="Firing rate [Hz]", colormap=colormap)
    lg1 = GridLayout(lg[2,5:6])        
    lscene6 = plot_response_fields!(lg1, rf;filter_spurious=true, colorbar_below=true, colormap=:rain,floor_offset=floor_offset, ceiling_offset=ceiling_offset)
    center!(lscene6.scene)
    fig = lg.parent.parent
    if do_animate
        display(fig)
        θ = 0
        fps = 10
        if fname === nothing
            mfname = "test4.mp4"
            do_wait = true 
        else
            b,e = splitext(fname)
            mfname = replace(fname, e => ".mp4")
            do_wait = false
        end
        record(fig.scene, mfname;framerate=fps,compression=10) do io
            while θ < 2π
                for lscene in [lscene1, lscene2, lscene3, lscene4, lscene5, lscene6]
                    rotate_cam!(lscene.scene, Vec3f(0.0, π/48,0.0))
                end
                θ += π/48
                if do_wait
                    sleep(1/fps)
                end
                recordframe!(io)
            end
        end
    else
        link_cameras_lscene(fig)
    end
end

function plot_spatial_responses(celldirs::Vector{String};kwargs...)
    h = process_kwargs(SpatialMapNew;kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    hs = string(h,base=16)
    outdir = joinpath(@__DIR__, "..","figures","single_cells", "spatial_response_$(hs)")
    if !isdir(outdir) 
        mkpath(outdir)
    end
    for celldir in celldirs
        cn = DPHT.get_shortname(celldir)
        fig = plot_spatial_responses(celldir;kwargs...)
        save(joinpath(outdir, "$(cn).png"),fig;px_per_unit=4)
    end
end


function plot_view_responses(celldirs::Vector{String};kwargs...)
    h = process_kwargs(ViewMapNew;kwargs...)
    h = process_kwargs(GazeResponseFields,h;kwargs...)
    hs = string(h,base=16)
    outdir = joinpath(@__DIR__, "..","figures","single_cells", "view_response_$(hs)")
    if !isdir(outdir) 
        mkpath(outdir)
    end
    for celldir in celldirs
        cn = DPHT.get_shortname(celldir)
        fname = joinpath(outdir, "$(cn).png") 
        fig = plot_view_responses(celldir;do_animate=true, fname=fname, kwargs...)
    end
end

function plot_spatial_maps(celldirs::Vector{String};use_response_field_rate=false, max_n_cols=10, _plot_theme=plot_theme, kwargs...)
    
    n = length(celldirs)
    nrows = round(Int64, ceil(max_n_cols))
    ncols = min(n, max_n_cols)
    @show nrows ncols
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=nrefinements.p));
    m_floor_hr = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));

    λ_all = map(celldirs) do celldir
        jm,rf = cd(celldir) do
            jm = JointMap(;kwargs...)
            rf = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 1000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001, use_trials=:all)
            jm, rf
        end
        spm = SpatialMapNew(jm,m_floor)
        if use_response_field_rate
            λ = rf.λ
        else
            λ=get_rate_map(spm)
        end
        # get an esimate of the size of the spatial fields
        nclusters = get_num_fields(rf)
        clusters = merge_fields(rf)
        cidx = findall(dropdims(mean(nclusters,dims=2),dims=2).< 0.001)
        spatial_fields = [rf.binidx[c] for c in clusters[cidx]]
        λ, spatial_fields
    end
    with_theme(_plot_theme) do 
        fig = Figure(size=(250*ncols, 250*nrows))
        for (i,(λ,ff)) in enumerate(λ_all)
            r = div(i,ncols)+1
            c = i - (r-1)*ncols
            ax = Axis(fig[r,c],aspect=1.0)
            hidedecorations!(ax)
            ax.bottomspinevisible = false
            ax.leftspinevisible = false
            viz!(ax, m_floor;color=:gray)
            if use_response_field_rate
                viz!(ax, m_floor_hr;color=λ, colormap=:rain)
            else
                viz!(ax, m_floor;color=λ, colormap=:rain)
            end
            for _ff in ff
                bb = find_boundary(m_floor_hr[_ff])
                viz!(ax, bb, color=:red)
            end
            if i == length(λ_all)
                Colorbar(fig[r,c+1], colorrange=(0,1), colormap=:rain, label="Firing rate", tellwidth=false, tellheight=false,alignmode=Outside(), ticksvisible = false, ticklabelsvisible=false)
            end

        end
        fig
    end
end

function wedge(r1, r2, θ1, θ2;origin=(0.0, 0.0))
    BezierPath([MoveTo(r1*cos(θ1), r1*sin(θ1)),
    EllipticalArc(origin[1], origin[2], r1, r1, 0.0, θ1, θ2),
    LineTo(r2*cos(θ2), r2*sin(θ2)),
    EllipticalArc(origin[1], origin[2], r2, r2, 0.0, θ2, θ1),
    LineTo(r1*cos(θ1), r1*sin(θ1)),
    ])
end

function plot_polar_histogram!(ax, r_data::AbstractVector{T}, θ_data::AbstractVector{T};nbins=(20,20)) where T <: Real
    # actually, we can just use regular hist for this
    θbins = range(-π, stop=π, length=nbins[1])
    rbins = range(0.0, stop=maximum(r_data), length=nbins[2])
    Z = zeros(length(θbins), length(rbins))
    for (r,θ) in zip(r_data, θ_data)
        i = searchsortedfirst(θbins, θ)
        j = searchsortedfirst(rbins, r)
        Z[i,j] += 1.0
    end
    plot_polar_histogram!(ax, Z, (θbins, rbins))
end

function plot_polar_histogram!(ax, Z::Matrix{<:Real}, bins::Tuple{T2, T2}) where T2 <: AbstractVector{T} where T <: Real
    # create wedges for the bins
    θbins,rbins = bins
    cr = extrema(Z)
    for (j,(r0,r1)) in enumerate(zip(rbins[1:end-1], rbins[2:end]))
        for (i,(θ0, θ1)) in enumerate(zip(θbins[1:end-1], θbins[2:end]))
            w = wedge(r0,r1, θ0, θ1)
            poly!(ax, w, color=Z[i,j], colorrange=cr)
        end
    end
end