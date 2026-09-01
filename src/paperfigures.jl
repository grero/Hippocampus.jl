module PaperFigures
using CairoMakie
using GLMakie
using Meshes
using StatsBase
using Hippocampus
using JLD2

plot_theme = Theme(Axis=(xlabelsize=14, ylabelsize=14,
                           xticklabelsize=14, yticklabelsize=14,
                           topspinevisible=false, rightspinevisible=false,
                           xgridvisible=false, ygridvisible=false,ylabelvisible=true,
                           xticklabelsvisible=true, xlabelvisible=true),
                     Scatter=(markersize=10px,),
                     Lines=(linewidth=3,),
                     fontsize=14)


function get_sessions()
    allcelldirs = open("/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/cell_list.txt") do fid
        readlines(fid)
    end
    sessiondirs = unique(Hippocampus.DPHT.get_level_path.("session", allcelldirs))
end

"""
    figure1()

Behavioural task and recording sites
"""
function figure1()
    udata = cd("/Volumes/Hippocampus/Data/picasso-misc/20181102/session01") do
        Hippocampus.UnityData()
    end
    poster_pos =  Dict(k=>v[[1,3]] for (k,v) in udata.header["PosterLocations"]) 
    sessions = get_sessions()
    with_theme(plot_theme) do
        fig = Figure(size=(768,768))
        lg_top = GridLayout(fig[1,1:2])
        axg = Axis(lg_top[1,1],aspect=DataAspect())
        Label(lg_top[1,1,TopLeft()], "A")
        hidedecorations!(axg)
        hidespines!(axg)
        img = load(joinpath(@__DIR__, "..","artefacts","monkey_in_chair.png"))
        image!(axg, rotr90(img))
        # load chamber rendering
        img_chamber = load(joinpath(@__DIR__, "..","figures","picasso_chamber_rendering.png"))
        axc = Axis(lg_top[1,2], aspect=DataAspect())
        Label(lg_top[1,2,TopLeft()],"E")
        hidedecorations!(axc)
        hidespines!(axc)
        image!(axc, rotr90(img_chamber))
        lg1 = GridLayout(fig[2,1])
        lg11 = GridLayout(lg1[1,1])
        Label(lg1[1,1,TopLeft()],"B")
        Hippocampus.plot_flat_maze_with_posters!(lg11;_poster_pos=poster_pos)
        lg12 = GridLayout(lg1[2,1])
        Label(lg1[2,1, TopLeft()], "C", padding=(20,20,0,0))
        trajectory_lengths = JLD2.load(joinpath(@__DIR__, "..", "data","optimal_path_length.jld2"), "lengths")
        Hippocampus.plot_trajectory_lengths!(lg12, trajectory_lengths, markersize=20)
        # show trajectories for a sample session
        lg2 = GridLayout(fig[2,2])
        Label(lg2[1,1,TopLeft()], "D")
        Hippocampus.plot_trajectories!(lg2, udata) 
        colsize!(fig.layout,1, Relative(0.3))

        fig
    end
end
"""
Place cells
"""
function figure2()
    spatially_selective, example_idx = JLD2.load(joinpath(@__DIR__, "..","data","paper","figure2_data.jld2"), "spatially_selective", "example_cell_idx")
    figure2(spatially_selective, example_idx)
end

function figure2(spatial_cells::Vector{String},example_idx::Vector{Int64})
    CairoMakie.activate!()
    plot_field_summary(Hippocampus.SpatialResponseFields, spatial_cells, example_idx)
end

function plot_response_field_with_sic(celldir::String,args...)
    with_theme(plot_theme) do
        fig = Figure(size=(600,600))
        lg = GridLayout(fig[1,1])
        plot_response_field_with_sic!(lg, celldir,args...)
        fig
    end
end

function plot_response_field_with_sic!(lg, celldir::String, ::Type{T}) where T <: Hippocampus.AbstractResponseFields
    kwargs = (min_speed=1.0, trial_start=2, min_place_obs=5, min_place_duration=0.05, min_view_obs=5, min_view_duration=0.02, pv_threshold=0.001)
    mm = Hippocampus.get_mesh(T,(p=3,g=2))
    sic,rfs = cd(celldir) do
        sic = Hippocampus.compute_skaggs_sic(Hippocampus.get_sic_type(T), 10_000;load_only=true, smooth=true, smoothing_method=:laplace, α=0.1, niter=50,kwargs...)
        rfs = Hippocampus.get_response_fields(T, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,pv_threshold=0.001,kwargs...) 
        sic,rfs
    end
    if isempty(rfs.binidx) || sic === nothing
        @show jj
    end
    lg12 = GridLayout(lg[1,1])
    Z = rfs.λ
    if embeddim(mm) == 2
        ax = Axis(lg12[1,1],aspect=1)
        hidedecorations!(ax)
        ax.leftspinevisible = false
        ax.bottomspinevisible = false
        # show the outline of the maze
        viz!(ax, mm;color=:lightgray)
        viz!(ax,mm;color=Z)
    else
        ax = LScene(lg12[1,1],show_axis=false)
        #TODO: only hide the ceiling if this cell has no field on the ceiling
        Hippocampus.plotmesh!(ax, mm;color=Z,indicate_north=true,hide_ceiling=true)
    end

    clusters = Hippocampus.merge_fields(rfs)
    nclusters = Hippocampus.get_num_fields(rfs)
    # relative number of times we get cluster of at least length.(clusters) randomly
    threshold = dropdims(sum(nclusters,dims=2),dims=2)/size(nclusters,2)
    # only keep fields where the probabilty of getting the same field in the surroages is less than 0.01
    valid_cluster_idx = findall(threshold .< kwargs.pv_threshold)
    for cluster in clusters[valid_cluster_idx] 
        bb = Hippocampus.find_boundary(mm, rfs.binidx[cluster])
        viz!(ax, bb;color=:black)
    end
    Colorbar(lg12[2,1], colorrange=extrema(filter(isfinite, Z)), label="Firing rate [Hz]",vertical=false,tellwidth=false)
    #colsize!(lg1, 1, Relative(0.6))

    # show SIC distribution
    #lg2 = GridLayout(fig[1,2])
    ax2 = Axis(lg12[3,1])

    boxplot!(ax2, fill(1.0, length(sic.sic)), sic.sic,color=:gray,show_notch=true,orientation=:horizontal,show_outliers=false)
    vlines!(ax2, sic.sic0, linestyle=:dot, color=:black)
    ax2.yticklabelsvisible = false
    ax2.yticksvisible = false
    ax2.bottomspinevisible = true 
    ax2.leftspinevisible = false
    ax2.rightspinevisible = false 
    ax2.xlabel = "SIC"
    rowsize!(lg12, 1, Relative(0.8))
end

function plot_field_summary(::Type{T}, celldirs::Vector{String},example_idx::Vector{Int64};plot_kwargs...) where T<:Hippocampus.AbstractResponseFields
    kwargs = (min_speed=1.0, trial_start=2, min_place_obs=5, min_place_duration=0.05, min_view_obs=5, min_view_duration=0.02, pv_threshold=0.001)
    mm = Hippocampus.get_mesh(T,(p=3,g=2))
    #m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3))
    with_theme(plot_theme) do
        fig = Figure(size=(900,800))
        # Example of a cell with place activity
        lgm = GridLayout(fig[1,1])
        Label(lgm[1,1,TopLeft()], "A")
        for (ii,jj) in enumerate(example_idx)
            lg1 = GridLayout(lgm[1,ii])
            sic,rfs = cd(celldirs[jj]) do
                sic = Hippocampus.compute_skaggs_sic(Hippocampus.get_sic_type(T), 10_000;load_only=true, smooth=true, smoothing_method=:laplace, α=0.1, niter=50,kwargs...)
                rfs = Hippocampus.get_response_fields(T, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,pv_threshold=0.001,kwargs...) 
                sic,rfs
            end
            if isempty(rfs.binidx) || sic === nothing
                @show jj
            end
            lg12 = GridLayout(lg1[1,1])
            show_points = T <: Hippocampus.GazeResponseFields
            show_boundaries = !show_points
            Hippocampus.plot_response_fields!(lg12, rfs;filter_spurious=true, colorbar_below=true, colormap=:jet, show_points=show_points, show_boundaries=show_boundaries, plot_kwargs...)

            # show SIC distribution
            #lg2 = GridLayout(fig[1,2])
            ax2 = Axis(lg12[3,1])
            if ii == 1
                Label(lg12[3,1, TopLeft()],"B";tellheight=false)
            else
                Label(lg12[3,1, TopLeft()],"")
            end
            boxplot!(ax2, fill(1.0, length(sic.sic)), sic.sic,color=:gray,show_notch=true,orientation=:horizontal,show_outliers=false)
            vlines!(ax2, sic.sic0, linestyle=:dot, color=:black)
            ax2.yticklabelsvisible = false
            ax2.yticksvisible = false
            ax2.bottomspinevisible = true 
            ax2.leftspinevisible = false
            ax2.rightspinevisible = false 
            ax2.xlabel = "SIC"
            rowsize!(lg12, 1, Relative(0.9))
        end

        # summary showing total number of fields and coverage
        lg3 = GridLayout(fig[2,1]) 
        Hippocampus.plot_n_fields!(lg3, T,celldirs;labels=["C","D","E","F"], pv_threshold=0.001, smooth=true, smoothing_method=:laplace, α=0.1,niter=50,colormap=:jet,floor_offset=-20, kwargs...)
        rowsize!(fig.layout, 1, Relative(0.6))
        fig
    end
end


"""
View cells
"""
function figure3(view_cells::Vector{String}, example_idx::Vector{<:Integer};kwargs...)
    # plot this using GLMakie
    GLMakie.activate!()
    plot_field_summary(Hippocampus.GazeResponseFields, view_cells, example_idx;hide_ceiling=true, indicate_north=false,show_points=true, show_boundaries=false, pointsize=2.5, floor_offset=-20.0, kwargs...)
end


"""
Mixed selective and conjunction cells
"""
function figure4()
end

function get_performance(session::String)
     udata = cd(session) do
        Hippocampus.UnityData()
    end
    perf = zeros(6)
    timeouts = zeros(6)
    for k in 1:6
        cidx = findall(udata.triggers[:,3].==30+k)
        # exclude repeat trials
        # that is trials preceded by an incorrect trial with the same poster
        cidx = filter(c->c==1 ? true : ((udata.triggers[c-1,3] .< 40)&&(udata.triggers[c-1,1]!=udata.triggers[c,1])), cidx)
        perf[k] = length(cidx)/sum(udata.triggers[:,1].==10+k)
        timeouts[k] = sum(udata.triggers[:,3].==40+k)/sum(udata.triggers[:,1].==10+k)
    end
    perf, timeouts
end
"""
Plot the performance of the animal
"""
function plot_performance()
    fname = joinpath(@__DIR__,"..", "data","performance.jld2")
    if isfile(fname)
        perf,timeouts = JLD2.load(fname, "performance","timeouts")
    else
        sessions = get_sessions()
        # for each session, get he performance
        # TODO: Filter repeat trials, i.e. trials in which the monkey failed and then repeated the same trial
        perf = zeros(6, length(sessions))
        timeouts = zeros(6, length(sessions))
        for (ii,session) in enumerate(sessions)
            perf[:,ii], timeouts[:,ii] = get_performance(session)
        end
        JLD2.save(fname, Dict("performance"=>perf, "timeouts"=>timeouts, "sessions"=>sessions))
    end
    @show mean(perf) std(perf) percentile(perf[:], [25, 50, 75])
    @show mean(timeouts) std(timeouts) percentile(timeouts[:], [10,50,95]) extrema(timeouts[:])
    imgs = [load(Hippocampus.poster_img[nn]) for nn in Hippocampus.poster_names]
    with_theme(plot_theme) do
        fig = Figure()
        ax = Axis(fig[1,1])
        ii = collect(CartesianIndices(size(perf)))
        boxplot!(ax, getindex.(ii, 1)[:], perf[:];show_outliers=false)
        ax.xticksvisible = false
        ax.xticklabelsvisible = false
        ax.bottomspinevisible = false
        # TODO: Use the posters as axis tick labels
        scatter!(ax, [1:6;], fill(0.65, 6), marker=imgs, markersize=70)
        ylims!(ax, 0.6, 1.0)
        ax.yticks = [0.7, 0.8, 0.9, 1.0]
        ax.ytrimspine = true
        fig
    end
end

function get_trajectory_time(udata::Hippocampus.UnityData)
    nt = Hippocampus.numtrials(udata)
    Δt = zeros(6, 6)
    Δt² = zeros(6, 6)
    nn = zeros(Int64, 6, 6)
    for i in 2:nt
        # make sure both the current and the previous trials were correct
        if  (30 .< udata.triggers[i,3] .< 40) && (30 .< udata.triggers[i-1,3] .< 40)
            k1 = udata.triggers[i-1,1] .- 10
            k2 = udata.triggers[i,1] .- 10
            _Δt = udata.timestamps[i,3] - udata.timestamps[i,2]
            Δt[k2,k1] += _Δt
            Δt²[k2,k1] += _Δt^2
            nn[k2,k1] += 1
        end
    end
    Δt ./= nn
    Δt² ./= nn
    Δt, Δt²
end

function get_trajectory_time(sessions::Vector{String};redo=false)
    fname = joinpath(@__DIR__, "..","data","trajectory_time.jld2")
    if !redo && isfile(fname)
        Δt,Δt² = JLD2.load(fname, "Δt","Δt²")
    else
        nn = length(sessions)
        Δt = zeros(6,6,nn)
        Δt² = zeros(6,6,nn)
        for (ii,session) in enumerate(sessions)
            udata = cd(session) do
                Hippocampus.UnityData()
            end
            Δt[:,:,ii],Δt²[:,:,ii] = get_trajectory_time(udata)
        end
        JLD2.save(fname, Dict("Δt"=>Δt, "Δt²"=>Δt², "sessions"=>sessions))
    end
    Δt, Δt²
end

function plot_trajectory_time!(lg, udata::Hippocampus.UnityData;kwargs...)
    Δt,Δt² = get_trajectory_time(udata) 
    plot_trajectory_time!(lg, sqrt(Δt².-Δt^2);kwargs...)
end

function plot_trajectory_time!(lg, Δt::Array{<:Real, 3};kwargs...)
    plot_trajectory_time!(lg, dropdims(mean(Δt, dims=3),dims=3);kwargs...)
end

function plot_trajectory_time!(lg, Δt::Matrix{<:Real};kwargs...)
    markersize = get(kwargs, :markersize, 45)
    ax = Axis(lg[1,1])
    h = heatmap!(ax, Δt)
    Colorbar(lg[1,2], h, label="CV(traj time)")

    imgs = [load(Hippocampus.poster_img[nn]) for nn in Hippocampus.poster_names]
    # create dummy axes for the labels
    axl = Axis(lg[1,0])
    scatter!(axl, fill(0.0, 6), [1:6;], marker=imgs, markersize=markersize)
    axb = Axis(lg[2,1])
    scatter!(axb, [1:6;], fill(0.0, 6), marker=imgs, markersize=markersize)
    colsize!(lg, 0, 2*markersize-20)
    rowsize!(lg, 2, 2*markersize-20)
    linkxaxes!(axb, ax)
    linkyaxes!(axl, ax)
    hidedecorations!(axb)
    hidespines!(axb)
    hidedecorations!(axl)
    hidespines!(axl)
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xticks = [1:6;]
    ax.yticks = [1:6;]
    axl.xlabel = "From"
    axl.xlabelvisible = true
    axb.ylabel = "To"
    axb.ylabelvisible = true
    colgap!(lg, 1, 1)
    rowgap!(lg, 1, 1)
end

function plot_trajectory_time(args...;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_trajectory_time!(lg, args...;kwargs...)
        fig
    end
end

function plot_directional_place_fields()     
    place_cells_width_directed_fields = readlines(open(joinpath(@__DIR__, "..","data/place_cells_with_oriented_fields.txt")))
    # place_cells_width_directed_fields[13]
    kwargs = (only_full_traversal=true, nshuffles=1000, nrefinements=(p=3,g=2), min_speed=1.0, trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, pv_threshold=0.001, redo=fname->false,colormap=:jet)
    celldir = "/Volumes/Hippocampus/Data/picasso-misc/20180727/session01/array03/channel086/cell01"
    with_theme(plot_theme) do
        fig = Figure(size=(700,550))
        lg1 = GridLayout(fig[1,1])
        lg2 = GridLayout(fig[2,1])
        rowsize!(fig.layout, 1, Relative(0.6))
        Hippocampus.plot_field_direction_tuning!(lg1, Hippocampus.MajorAxisDirectionTuning, celldir, nothing;start_label='A', kwargs...)
        Hippocampus.plot_field_direction_summary!(lg2, Hippocampus.MajorAxisDirectionTuning, place_cells_width_directed_fields) 
        Label(lg2[1,1, TopLeft()], "E")
        fig
    end
end

function plot_simpler_place_cells(;figsize=(600,600),_plot_theme=plot_theme)
    ff,args = JLD2.load(joinpath(@__DIR__, "..", "data", "simpler_place_field_estimation.jld2"), "field_index","kwargs")
    nrefinements = args["nrefinements"]
    n_fields = length.(ff)
    m_floor = Hippocampus.floor_topology3(;nrefinements=nrefinements.p)
    # compute coverage
    Z = zeros(nelements(m_floor))
    for _ff in ff
        for __ff in _ff
            Z[__ff] .+= 1.0
        end
    end

    with_theme(_plot_theme) do
        fig = Figure(size=figsize)
        ax1 = Axis(fig[1,1])

        hist!(ax1, n_fields, color=:gray45)
        Label(fig[1,1,TopLeft()], "A")
        ax1.xticks = n_fields
        ax1.xlabel = "Number of fields"
        ax1.ylabel = "Number of cells"

        #coverage
        ax2 = Axis(fig[1,2], aspect=1)
        Label(fig[1,2,TopLeft()], "B")
        hidedecorations!(ax2)
        viz!(ax2, m_floor;color=:lightgray)
        viz!(ax2, m_floor;color=Z, colormap=:jet)
        hidespines!(ax2)
        Hippocampus.plot_pillars!(ax2)
        fig
    end
end

end #module