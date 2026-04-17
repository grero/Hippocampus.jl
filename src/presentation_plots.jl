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

function plot_directional_conditioned_view_field_example(celldir::String;kwargs...)
     rf, rfsp, gidx = cd(celldir) do 
        rf = Hippocampus.get_response_fields(Hippocampus.GazeResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        rfsp = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        gidx = Hippocampus.DirectionFiltered(;only_full_traversal=true, min_speed=1.0, trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, pv_threshold=0.001,redo=fname->false)
        rf,rfsp, gidx
    end
    μr,ϕ = Hippocampus.get_directional_tuning_strength(gidx;do_shuffle=false, smooth=false,niter=1)
    res,vc,λ = Hippocampus.analyse_directionality(gidx, rf;pv_threshold=0.01,niter=50)
    spatial_clusters = Hippocampus.merge_fields(rfsp)
    m_floor = Translate(0.0, 0.0, -20)(Hippocampus.floor_topology3(;nrefinements=3))
    cm = Point3f.(Tuple.(centroid.(m_floor[rfsp.binidx[spatial_clusters[1]]])))
    cmp = mean(cm)
    with_theme(_plot_theme)  do
        fig = Figure(size=(600,300))
        lg1 = GridLayout(fig[1,1])
        Hippocampus.plot_response_fields!(lg1, rf;colormap=:rain, filter_spurious=true, mazecolor=:darkgray,floor_offset=-20, ceiling_offset=15)
        lg2 = GridLayout(fig[1,2])
        lscene = Hippocampus.plot_response_fields!(lg2, rf, λ[1];colormap=:rain, filter_spurious=true, mazecolor=:darkgray, floor_offset=-20, ceiling_offset=15)
        #viz!(lscene, cm, color=:orange)
        scatter!(lscene, cm, color=:orange, markersize=10px)
        arrows3d!(lscene, cmp + Point3f(0.0, 0.0, 0.1), Point3f(2*cos(ϕ[1]),2*sin(ϕ[1]), 0.0),color=:black)
        Hippocampus.link_cameras_lscene(fig)
        fig
    end
end

function plot_view_conditioned_on_place(celldir::String;kwargs...)
     rf_gaze, rf_spatial,jm = cd(celldir) do 
        rf = Hippocampus.get_response_fields(Hippocampus.GazeResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        rfsp = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        jm = Hippocampus.JointMap(;redo=fname->true, do_save=false, nrefinements=(p=3,g=2),min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, trial_start=2)
        rf,rfsp, jm
    end
    view_clusters = Hippocampus.merge_fields(rf_gaze)
    spatial_clusters = Hippocampus.merge_fields(rf_spatial)
    
    res = Hippocampus.conjunctions(jm, rf_spatial,rf_gaze.binidx[view_clusters[1]])
    fidx = findall(isfinite, res[1]["λ_sub"])
    qidx = fidx[argmax(res[1]["λ_sub"][fidx])]

    vm = Hippocampus.ViewMapNew(jm, mm;placebins=res[1]["sub_idx"][qidx])
    vm_orig = Hippocampus.ViewMapNew(jm, mm;placebins=rf_spatial.binidx[spatial_clusters[1]]);
    m_floor = Translate(0.0, 0.0, -20)(Hippocampus.floor_topology3(;nrefinements=3))
    with_theme(_plot_theme) do
        fig = Figure(size=(600,500))
        lg1 = GridLayout(fig[1,1])
        lscene1 = Hippocampus.plot_response_fields!(lg1, rf_gaze, Hippocampus.get_rate_map(vm_orig), colormap=:rain, filter_spurious=true,mazecolor=:darkgray,floor_offset=-20, ceiling_offset=15)
        scatter!(lscene1, Point3f.(Tuple.(centroid.(m_floor[rf_spatial.binidx[spatial_clusters[1]]]))), color=:orange)
        lg2 = GridLayout(fig[1,2])
        lscene2 = Hippocampus.plot_response_fields!(lg2, rf_gaze, Hippocampus.get_rate_map(vm), colormap=:rain, filter_spurious=true, mazecolor=:darkgray, floor_offset=-20, ceiling_offset=15)
        scatter!(lscene2, Point3f.(Tuple.(centroid.(m_floor[res[1]["sub_idx"][qidx]]))), color=:orange)
        Hippocampus.link_cameras_lscene(fig)
        fig
    end

end

function plot_cell_categories(;kwargs...)
    allcelldirs = open("/Volumes/Hippocampus/Data/picasso-misc/AnalysisHM/Current Analysis/cell_list.txt") do fid
       readlines(fid)
    end
    spatially_selective = JLD2.load(joinpath(@__DIR__, "..","data","paper","figure2_data.jld2"), "spatially_selective")
    view_selective = JLD2.load(joinpath(@__DIR__, "..","data","paper","figure3_data.jld2"), "view_selective")
    non_selective = setdiff(allcelldirs, union(spatially_selective, view_selective))
    num_place_fields = JLD2.load(joinpath(@__DIR__, "..","field_stats_73eba60a.jld2"), "nfields")
    num_view_fields = JLD2.load(joinpath(@__DIR__, "..","field_stats_a91f3809.jld2"), "nfields")
    place_cells = collect(keys(filter(k->k[2]>0, num_place_fields)))
    directional_tuning_results = JLD2.load(joinpath(@__DIR__,"..","data","place_field_directional_cells.jld2"),"directional_tuning_results")
    directional_place_cells = place_cells[directional_tuning_results.==true]
    view_cells = collect(keys(filter(k->k[2]>0, num_view_fields)))
    place_and_view_selective = intersect(spatially_selective, view_selective)
    res_pvc = JLD2.load("data/place_accounting_for_view.jld2", "pvc_results")
    conjunctive = place_and_view_selective[res_pvc.==true]
    colors = Dict(:directional => :yellow, :place_cells => :orange, spatially_selective => :red,
                  :view_selective => :blue, :view_cells => :green, :conjunctive => :purple)

    with_theme(_plot_theme) do
        fig = Figure(size=(300,300))
        ax = Axis(fig[1,1])
        barplot!(ax, [1:8;], length.([non_selective, spatially_selective, place_cells, directional_place_cells, view_selective, view_cells, place_and_view_selective, conjunctive]), 
                                    color=[:gray, :red, :orange, :yellow, :blue, :green,:purple, :purple4])
        ax.bottomspinevisible = false 
        ax.xticklabelsvisible = false
        ax.ylabel = "Number of cells"
        fig
    end
end