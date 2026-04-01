using CairoMakie
using Hippocampus

_plot_theme = theme_dark()
_plot_theme.Axis.xgridvisible = false
_plot_theme.Axis.ygridvisible = false
_plot_theme.Axis.leftspinevisible = true
_plot_theme.Axis.leftspinecolor = :gray45
_plot_theme.Axis.yticksvisible = true
_plot_theme.Axis.ytickcolor = :gray45
_plot_theme.Colorbar.ticksvisible = true
_plot_theme.Colorbar.tickcolor = :gray45

function plot_place_field_example(celldir::String;kwargs...)
    rf_spatial,sic = cd(celldir) do
        rf = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        sic =  Hippocampus.compute_skaggs_sic(Hippocampus.SpatialInformationContent, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,load_only=false,redo=fname->false,trial_start=2, min_speed=1.0, min_place_obs=5, min_place_duration=0.05, min_view_duration=0.02, min_view_obs=5)
        rf,sic
    end
    with_theme(_plot_theme) do
        fig = Figure(size=(600,350))
        lg1 = GridLayout(fig[1,1])
        Hippocampus.plot_response_fields!(lg1, rf_spatial;colormap=:rain,mazecolor=:darkgray)
        # plot sic
        ax = Axis(fig[1,2])
        boxplot!(ax, fill(1.0, length(sic.sic)),sic.sic, show_outliers=false, show_notch=true, color=:gray)
        hlines!(ax, sic.sic0, linestyle=:dot, color=:gray45)
        ax.bottomspinevisible = false
        ax.xticklabelsvisible = false
        ax.xticksvisible = false
        ax.ylabel = "SIC"
        colsize!(fig.layout, 2, Relative(0.2))
        fig
    end
end

function plot_view_field_example(celldir::String;kwargs...)
    rf_gaze,sic = cd(celldir) do
        rf = Hippocampus.get_response_fields(Hippocampus.GazeResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        sic =  Hippocampus.compute_skaggs_sic(Hippocampus.GazeInformationContent, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,load_only=false,redo=fname->false,trial_start=2, min_speed=1.0, min_place_obs=5, min_place_duration=0.05, min_view_duration=0.02, min_view_obs=5)
        rf,sic
    end
    with_theme(_plot_theme) do
        fig = Figure(size=(600,350))
        lg1 = GridLayout(fig[1,1])
        Hippocampus.plot_response_fields!(lg1, rf_gaze;colormap=:rain,kwargs...)
        # plot sic
        ax = Axis(fig[1,2])
        boxplot!(ax, fill(1.0, length(sic.sic)),sic.sic, show_outliers=false, show_notch=true, color=:gray)
        hlines!(ax, sic.sic0, linestyle=:dot, color=:gray45)
        ax.bottomspinevisible = false
        ax.xticklabelsvisible = false
        ax.xticksvisible = false
        ax.ylabel = "SIC"
        colsize!(fig.layout, 2, Relative(0.2))
        fig
    end
end


function plot_directionality_example(celldir::String;kwargs...)
    rf, gidx = cd(celldir) do 
        rf = Hippocampus.get_response_fields(Hippocampus.GazeResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        gidx = Hippocampus.DirectionFiltered(;only_full_traversal=true, min_speed=1.0, trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, pv_threshold=0.001,redo=fname->false)
        rf,gidx
    end
    # get the behavioural data as well
    sessiondir = Hippocampus.DPHT.get_level_path("session", celldir)
    udata = cd(sessiondir) do
        Hippocampus.UnityData()
    end
    with_theme(_plot_theme) do
        fig = Figure(size=(600,300))
        lg1 = GridLayout(fig[1,1])
        Hippocampus.plot_response_fields!(lg1, rf_spatial;colormap=:rain,mazecolor=:darkgray,colorbar_below=true)
        lg2 = GridLayout(fig[1,2])
        lg21 = GridLayout(lg2[1,1])
        Hippocampus.plot_directional_tuning!(lg21, gidx;kwargs...)
        Label(lg2[1,2], "Firing rate", rotation=π/2, tellheight=false)
        lg22 = GridLayout(lg2[2,1])
        Hippocampus.plot_directional_tuning!(lg22, gidx;occupancy_only=true, kwargs...)
        Label(lg2[2,2], "Occupancy", rotation=π/2, tellheight=false)
        lg3 = GridLayout(fig[1,3])
        Hippocampus.plot_field_traversals!(lg3, gidx, udata;floorcolor=:darkgray)
        colsize!(fig.layout, 2, Relative(0.2))
        fig
    end
end