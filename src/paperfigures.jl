module PaperFigures
using CairoMakie
using GLMakie
using Meshes
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

"""
    figure1()

Behavioural task and recording sites
"""
function figure1()
    udata = cd("/Volumes/Hippocampus/Data/picasso-misc/20181102/session01") do
        Hippocampus.UnityData()
    end
    with_theme(plot_theme) do
        fig = Figure(size=(768,358))
        lg1 = GridLayout(fig[1,1])
        Hippocampus.plot_flat_maze_with_posters!(lg1)
        # show trajectories for a sample session
        lg2 = GridLayout(fig[1,2])
        Hippocampus.plot_trajectories!(lg2, udata) 
        colsize!(fig.layout,1, Relative(0.35))
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

end #module