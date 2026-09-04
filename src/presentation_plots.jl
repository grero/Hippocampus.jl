using CairoMakie
using Hippocampus
using Random
using Makie.Colors

_plot_theme = theme_dark()
_plot_theme.Axis.xgridvisible = false
_plot_theme.Axis.ygridvisible = false
_plot_theme.Axis.leftspinevisible = true
_plot_theme.Axis.bottomspinevisible = true
_plot_theme.Axis.leftspinecolor = :gray45
_plot_theme.Axis.bottomspinecolor = :gray45
_plot_theme.Axis.yticksvisible = true
_plot_theme.Axis.xticksvisible = true
_plot_theme.Axis.ytickcolor = :gray45
_plot_theme.Axis.xtickcolor = :gray45
_plot_theme.Colorbar.ticksvisible = true
_plot_theme.Colorbar.tickcolor = :gray45
_plot_theme.fontsize = 18

function plot_place_field_example(celldir::String;kwargs...)
    rf_spatial,sic = cd(celldir) do
        rf = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0,trial_start=2, pv_threshold=0.001)
        sic =  Hippocampus.compute_skaggs_sic(Hippocampus.SpatialInformationContent, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,load_only=false,redo=fname->false,trial_start=2, min_speed=1.0, min_place_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, min_view_obs=-1)
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
        rf = Hippocampus.get_response_fields(Hippocampus.GazeResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0,trial_start=2, pv_threshold=0.001)
        sic =  Hippocampus.compute_skaggs_sic(Hippocampus.GazeInformationContent, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=50,load_only=false,redo=fname->false,trial_start=2, min_speed=1.0, min_place_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, min_view_obs=-1)
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

function plot_direction_selective_model()
    # model direcion selectivity
    X,X1,X2,Y, λ, northward_trajectory, southward_trajectory = Hippocampus.simulate_fake_directionality(;additive=2.0, multiplicative=2.0) 
    # apply mild smoothing to avoid singularities
    Xs = permutedims(Hippocampus.laplace_smoothing(permutedims(X), Ls, 0.1;niter=50));
    X1s = permutedims(Hippocampus.laplace_smoothing(permutedims(X1), Ls, 0.1;niter=50));
    X2s = permutedims(Hippocampus.laplace_smoothing(permutedims(X2), Ls, 0.1;niter=50))
    Ys = permutedims(Hippocampus.laplace_smoothing(permutedims(Y), Ls, 0.1;niter=50));

     m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3));
    with_theme(_plot_theme) do
        fig = Figure(size=(1000,600))
        # show the real field, with directions superimposed
        ax1 = Axis(fig[1,1], aspect=1)
        Label(fig[1,1,TopLeft()], "A")
        viz!(ax1, m_floor, color=λ, colormap=:thermal)
        hidedecorations!(ax1)
        ax1.bottomspinevisible = false
        ax1.leftspinevisible = false
        # superimpose trajectories
        viz!(ax1, centroid.(m_floor[northward_trajectory]), color=:red)
        viz!(ax1, centroid.(m_floor[southward_trajectory]), color=:orange)
        arrows2d!(ax1, Point2f(0.0), Point2f(0.0, 1.0), color=:red)
        arrows2d!(ax1, Point2f(2.0, 1.0), Point2f(0.0, -1.0), color=:orange)
        lg = GridLayout(fig[1,2])
        Label(lg[1,1,TopLeft()], "B")
        ax2 = Axis(lg[1,1])
        lines!(ax2, 1:length(northward_trajectory), λ[northward_trajectory], color=:red)
        lines!(ax2, 1:length(southward_trajectory), λ[southward_trajectory], color=:orange)
        ax2.ylabel = "Firing rate [Hz]"
        # additive
        ax3 = Axis(lg[1,2])
        lines!(ax3, 1:length(northward_trajectory), 0.08*λ[northward_trajectory] .+ 2.0, color=:red)
        lines!(ax3, 1:length(southward_trajectory), λ[southward_trajectory], color=:orange)

        ax4 = Axis(lg[1,3])
        # multiplicative
        lines!(ax4, 1:length(northward_trajectory), 2.0*λ[northward_trajectory], color=:red)
        lines!(ax4, 1:length(southward_trajectory), λ[southward_trajectory], color=:orange)

        Label(lg[2,1,TopLeft()],"C")
        ax5 = Axis(lg[2,1])
        barplot!(ax5, [1:2;], dropdims(mean(Xs,dims=1),dims=1), color=[:red, :orange])
        ax5.ylabel = "Mean spike count"
        ax5.xticks = ([1,2], ["Northward", "Southward"])

        ax6 = Axis(lg[2,2])
        barplot!(ax6, [1:2;], dropdims(mean(X1s,dims=1),dims=1), color=[:red, :orange])
        ax6.xticks = ([1,2], ["Northward", "Southward"])

        ax7 = Axis(lg[2,3])
        barplot!(ax7, [1:2;], dropdims(mean(X2s,dims=1),dims=1), color=[:red, :orange])
        ax7.xticks = ([1,2], ["Northward", "Southward"])
        linkaxes!(ax5,ax6,ax7)
        for ax in [ax5,ax6,ax7]
            ax.xticklabelrotation = -π/6
        end
        # show joint information
        ee1 = Hippocampus.get_conditional_information(Xs,Ys)
        ee2 = Hippocampus.get_conditional_information(X1s,Ys)
        ee3 = Hippocampus.get_conditional_information(X2s,Ys)
        ax8 = Axis(lg[3,1:3])
        Label(lg[3,1, TopLeft()], "D")
        lines!(ax8, [1:3;], [ee1, ee2, ee3])
        scatter!(ax8, [1:3;], [ee1, ee2, ee3])
        ax8.xticks = ([1:3;], ["Place","Add","Mult"])
        ax8.ylabel = "I(S;D|P)"
        colsize!(fig.layout,1,Relative(0.6))
        fig
    end
end


function plot_direction_selectivity_example()
    datadir = "/Volumes/Hippocampus/Data/picasso-misc"
    celldir = joinpath(datadir ,"20180802","session01","array02","channel045","cell01")
    card = cd(celldir) do
       Hippocampus.CardinalPlaceFieldDirectionality(;nshuffles=1000,only_full_traversal=true, nrefinements=(p=3,g=2), min_speed=1.0, trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, pv_threshold=0.001,redo=fname->false)
    end
    @show Hippocampus.issignificant(card)
    @show Hippocampus.DPHT.get_shortname(celldir)
    fig = Hippocampus.plot_directional_place_map(card,1;_plot_theme=_plot_theme)

end

function plot_direction_selectivity_summary()
    larger, smaller, both = JLD2.load(joinpath(@__DIR__, "..","data","directional_selectivity_summary.jld2"),"larger","lower","both")
    none = (~).(larger.|smaller.|both)
    with_theme(_plot_theme) do
        fig = Figure(size=(400,400))
        ax = Axis(fig[1,1])
        barplot!(ax, 1:4, sum.([none, both, larger, smaller]))
        ax.xticks = (1:4, ["Non-\nselective","Both","larger","smaller"])
        ax.ylabel = "Number of cells"
        fig
    end
end

function plot_place_view_summary()
    larger, smaller, both = JLD2.load(joinpath(@__DIR__, "..","data","place_conditioned_on_view_summary.jld2"),"larger","lower","both")
    none = (~).(larger.|smaller.|both)
    with_theme(_plot_theme) do
        fig = Figure(size=(400,400))
        ax = Axis(fig[1,1])
        barplot!(ax, 1:4, sum.([none, both, larger, smaller]))
        ax.xticks = (1:4, ["Non-\nselective","Both","larger","smaller"])
        ax.ylabel = "Number of cells"
        fig
    end

end

function plot_model_cell_directional_place_field()
    testdata_dir = joinpath(@__DIR__,"..","data","testdata","ModelSubject","20260422")
    test_sessiondir = joinpath(testdata_dir, "session01")
    testcell = joinpath(sessiondir, "array01/channel001/cell03")
    rf_spatial,gidx,vpvrp = cd(testcell) do
        rf = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;nrefinements=(p=3,g=2),smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0,trial_start=2, pv_threshold=0.001)
        gidx = Hippocampus.DirectionFiltered(;only_full_traversal=true, nrefinements=(p=3,g=2), min_speed=1.0, trial_start=2, smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, pv_threshold=0.001,redo=fname->false, do_save=true)
        vpvrp = Hippocampus.ViewAndPlaceRepresentationNew(;redo=fname->false,do_save=true,trial_start=2)
        rf,gidx,vpvrp
    end
    jocc,qdata = cd(test_sessiondir)  do
        jocc = Hippocampus.JointOccupancy(;redo=false, do_save=false, nrefinements=(p=3, g=2),trial_start=2,min_speed=1.0)
        qdata = Hippocampus.UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
        jocc, qdata
    end
    (spike_count_1, occupancy_1), (spike_count_2, occupancy_2) = Hippocampus.get_major_axis_direction_tuning(gidx, rf_spatial, jocc, vpvrp,qdata)
    λ1 = spike_count_1./occupancy_1
    λ2 = spike_count_2./occupancy_2
    q1,q2 = Hippocampus.permutation_test(mean, λ1, λ2)
    v = Hippocampus.get_major_axis(rf_spatial)
    θ = atan(v[2,1], v[1,1])
    bidx1 = findall(cos.(gidx.anglebins .- θ) .> cos(π/6))
    trialidx1 = unique(getindex.(filter(q->in(bidx1)(q[4]), gidx.index[1]),5))
    bidx2 = findall(cos.(gidx.anglebins .- (θ-π)) .> cos(π/6))
    trialidx2 = unique(getindex.(filter(q->in(bidx2)(q[4]), gidx.index[1]),5))
    jm1, jm2 = cd(testcell) do
        jm1 = Hippocampus.JointMap(;redo=fname->false, do_save=false, nrefinements=(p=3,g=2),min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, trial_start=2, use_trials=trialidx1)
        jm2 = Hippocampus.JointMap(;redo=fname->false, do_save=false, nrefinements=(p=3,g=2),min_speed=1.0, min_place_obs=-1, min_view_obs=-1, min_place_duration=-1.0, min_view_duration=-1.0, trial_start=2, use_trials=trialidx2)
        jm1, jm2
    end
    spm1 = Hippocampus.SpatialMapNew(jm1,m_floor);
    spml_1 = Hippocampus.SmoothedMap(spm1;method=:laplace, α=0.1, niter=50)
    spm2 = Hippocampus.SpatialMapNew(jm2,m_floor);
    spml_2 = Hippocampus.SmoothedMap(spm2;method=:laplace, α=0.1, niter=50)
    with_theme(_plot_theme) do
        fig = Figure(size=(1100,400))
        lg1 = GridLayout(fig[1,1])
        ax = Hippocampus.plot_response_fields!(lg1, rf_spatial;_plot_theme=_plot_theme,colormap=:rain)
        arrows2d!(ax, Point2f(0.0, 0.0), Point2f(2.0*cos(θ), 2.0*sin(θ)), color=:black)
        lg2 = GridLayout(fig[1,2])
        Label(lg2[1,1,Top()], "Forward")
        ax2 = Hippocampus.plot_response_fields!(lg2, rf_spatial,Hippocampus.get_rate_map(spml_1);_plot_theme=_plot_theme,colormap=:rain, show_points=false)
        lg3 = GridLayout(fig[1,3])
        Label(lg3[1,1,Top()], "Backward")
        ax3 = Hippocampus.plot_response_fields!(lg3, rf_spatial,Hippocampus.get_rate_map(spml_2);_plot_theme=_plot_theme,colormap=:rain, show_points=false)
        lg4 = GridLayout(fig[1,4])
        ax4 = Axis(lg4[1,1])
        boxplot!(ax4, fill(1.0, length(q1)), q1-q2)
        hlines!(ax4, mean(λ1) - mean(λ2), color=:gray45, linestyle=:dot)
        ax4.bottomspinevisible = false
        ax4.xticklabelsvisible = false
        ax4.ylabel = "λ_forward - λ_backward"
        colsize!(fig.layout, 4, Relative(0.1))
        fig
    end
end

"""
    sample_map(pvc::Hippocampus.PlaceViewConjunction, vidx::Integer, pidx::Integer;nruns=200)

Sample random points from the view conditioned map in `pvc.λ_infield[:, vidx]`. 
"""
function sample_map(pvc::Hippocampus.PlaceViewConjunction, vidx::Integer, pidx::Integer;nruns=200)
    fidx = findall(isfinite.(pvc.λ_infield[:,vidx]))
    qidx = Hippocampus.getfields(pvc.spatial_fields;cluster_threshold=0.001)[pidx]
    μ = mean(pvc.spatial_fields.λ[qidx])
    μc = mean(filter(isfinite, pvc.λ_infield[qidx,vidx]))
    widx = setdiff(fidx, qidx)
    X = zeros(nruns)
    for i in 1:length(X)
        _idx = shuffle(widx)[1:length(qidx)]
        X[i] = mean(pvc.λ_infield[_idx,vidx])
    end
    X, μ, μc
end

function sample_map(pvc::Hippocampus.ViewPlaceConjunction, vidx::Integer, pidx::Integer;nruns=200)
    fidx = findall(isfinite.(pvc.λ_infield[:,pidx]))
    qidx = Hippocampus.getfields(pvc.view_fields; cluster_threshold=0.001)[vidx]

    μ = mean(pvc.view_fields.λ[qidx])
    μc = mean(filter(isfinite, pvc.λ_infield[qidx,pidx]))
    widx = setdiff(fidx, qidx)
    X = zeros(nruns)
    for i in 1:length(X)
        _idx = shuffle(widx)[1:length(qidx)]
        X[i] = mean(pvc.λ_infield[_idx,pidx])
    end
    X, μ, μc
end

function plot_conjunctions_new(celldir::String,args...;kwargs...)
    plot_conjunctions_new(Hippocampus.SpatialResponseFields, Hippocampus.GazeResponseFields, celldir, args...;kwargs...)
end

function plot_conjunctions_new(::Type{T1}, ::Type{T2}, celldir::String,spatial_field_idx::Integer, view_field_idx::Integer;_plot_theme=_plot_theme, distr_plot=:boxplot, kwargs...) where T1 <: Hippocampus.SpatialResponseFieldsAll where T2 <: Hippocampus.GazeResponseFieldsAll
    jm,vpc, pvc = cd(celldir) do
        jm = Hippocampus.JointMap(;nrefinements=(p=3,g=2),min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, trial_start=2)
        vpc = Hippocampus.ViewPlaceConjunction(T2,T1;nshuffles=1000, nrefinements=(p=3,g=2), smooth=true, smoothing_method=:laplace, α=0.1, niter=50, min_view_obs=5, min_place_obs=5, min_view_duration=0.02, min_place_duration=0.05, min_speed=1.0, pv_threshold=0.001)
        pvc = Hippocampus.PlaceViewConjunction(T1,T2;nshuffles=1000, nrefinements=(p=3, g=2), smooth=true, smoothing_method=:laplace, α=0.1, niter=50,  min_view_obs=5, min_place_obs=5, min_view_duration=0.02,  min_place_duration=0.05, min_speed=1.0, pv_threshold=0.001)
        jm, vpc, pvc
     end
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    m_floor = Hippocampus.repair_mesh(Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3)))
    pidx = Hippocampus.getfields(vpc.spatial_fields;cluster_threshold=0.001)[spatial_field_idx]
    vidx = Hippocampus.getfields(vpc.view_fields;cluster_threshold=0.001)[view_field_idx]
    vm = Hippocampus.get_view_rate_map(jm, pidx, size(vpc.λ_infield,1))
    λ_infield_view = vec(Hippocampus.laplace_smoothing(vm, mm, 0.1;niter=50))
    λ_infield_view[vm.occupancy.==0] .= NaN
    λ_infield_place = vec(Hippocampus.laplace_smoothing(Hippocampus.get_spatial_rate_map(jm, vidx, size(pvc.λ_infield,1)),m_floor,0.1;niter=50))
    X_place_nc, μ_place, μ_place_view = sample_map(pvc, view_field_idx, spatial_field_idx)
    X_view_nc, μ_view, μ_view_place = sample_map(vpc, view_field_idx,spatial_field_idx)
    colormap = get(kwargs, :colormap, :rain)
    ccolors = Hippocampus.get_colors(colormap)

    with_theme(_plot_theme) do
        fig = Figure(size=(1000,550))
        # show both place and view maps
        lg1 = GridLayout(fig[1,1])
        Label(lg1[1,1,TopLeft()], "A")
        lg1_1 = GridLayout(lg1[1,1])
        # view field
        Hippocampus.plot_response_fields!(lg1_1, pvc.view_fields;colormap=colormap, hide_ceiling=true, indicate_north=false, floor_offset=0.0)
        # place field
        lg1_2 = GridLayout(lg1[2,1])
        Label(lg1[2,1,TopLeft()],"B")
        Hippocampus.plot_response_fields!(lg1_2, pvc.spatial_fields;colormap=colormap)
        # view conditioned on place
        lg2 = GridLayout(fig[1,2])
        Label(lg2[1,1,TopLeft()], "C")
        lg2_1 = GridLayout(lg2[1,1])
        lscene = LScene(lg2_1[1,1], show_axis=false)
        Hippocampus.plot_pillars!(lscene)
        α_view = fill(0.0, length(λ_infield_view))
        α_view[isfinite.(λ_infield_view)] .= 1.0
        Hippocampus.plotmesh!(lscene, mm;color=λ_infield_view,alpha=α_view, colormap=colormap, indicate_north=false, hide_ceiling=true, showsegments=true,floor_offset=0.0)
        ax1 = Axis(lg2_1[2,1], aspect=1)
        hidedecorations!(ax1)
        Hippocampus.plot_pillars!(ax1)
        Z_floor = fill(parse(Colorant, :lightgray), nelements(m_floor))
        alpha = fill(0.0, nelements(m_floor))
        Z_floor[pidx] .= parse(Colorant, ccolors[spatial_field_idx])
        alpha[pidx] .= 1.0
        Z_floor[vpc.non_covered_idx[spatial_field_idx, view_field_idx]] .= parse(Colorant, :gray25)
        alpha[vpc.non_covered_idx[spatial_field_idx, view_field_idx]] .= 1.0
        alpha = fill(1.0, size(Z_floor)...)
        viz!(ax1, m_floor; color=Z_floor,alpha=alpha)
        # distribution of conditional and joint firing rate
        ax2 = Axis(lg2_1[3,1])
        Label(lg2_1[3,1, TopLeft()], "D")
        xx = [fill(1.0, length(X_view_nc));fill(2.0, length(vpc.λ_sub[spatial_field_idx,view_field_idx,:]))]
        yy = [X_view_nc;vpc.λ_sub[spatial_field_idx,view_field_idx,:]]
        ss_vpc = vpc.λ_covered[spatial_field_idx, view_field_idx] > percentile(vpc.λ_sub[spatial_field_idx, view_field_idx,:], 97.5)
        ss_vp = μ_view_place > percentile(X_view_nc, 97.5)
        cc = [fill(1, length(X_view_nc));fill(2, length(vpc.λ_sub[spatial_field_idx,view_field_idx,:]))]

        if distr_plot == :boxplot
            boxplot!(ax2, xx,yy, show_outliers=false,orientation=:horizontal, color=cc, colormap=[:royalblue, :gray25],show_notch=true)
        else
            violin!(ax2, xx,yy, orientation=:horizontal, show_median=true, color=[:royalblue, :gray25][cc])
        end
        if ss_vpc
            marker_vpc = :star6
        else
            marker_vpc = :circle
        end
        if ss_vp
            marker_vp = :start6
        else
            marker_vp = :circle
        end
        scatter!(ax2, [μ_view_place], [1.0],color=ccolors[view_field_idx], marker=marker_vp, label="Inside VF")
        scatter!(ax2, [μ_view], [1.0],color=:seagreen, label="Original VF")
        scatter!(ax2, [vpc.λ_covered[spatial_field_idx,view_field_idx]], [2.0], color=ccolors[spatial_field_idx], marker=marker_vpc, label="VF & PF")
        if ss_vpc
        end
        rowsize!(lg2_1, 3, 75)
        rowsize!(lg2_1,1, Relative(0.5))
        ax2.yticklabelsvisible = false
        ax2.yticksvisible = false
        ax2.xlabel = "Firing rate [Hz]"

        # place conditioned on view
        lg3 = GridLayout(fig[1,3])
        Label(lg3[1,1,TopLeft()], "E")
        lscene = LScene(lg3[1,1], show_axis=false)
        Hippocampus.plot_pillars!(lscene)
        Z_maze = fill(parse(Colorant, :lightgray), nelements(mm))
        alpha = fill(0.0, nelements(mm))
        Z_maze[vidx] .= parse(Colorant, ccolors[view_field_idx])
        alpha[vidx] .= 1.0
        Z_maze[pvc.non_covered_idx[spatial_field_idx, view_field_idx]] .= parse(Colorant, :gray25)
        alpha[pvc.non_covered_idx[spatial_field_idx, view_field_idx]] .= 1.0
        Hippocampus.plotmesh!(lscene, mm, color=Z_maze, alpha=alpha, indicate_north=false,hide_ceiling=true, showsegments=true)
        ax3 = Axis(lg3[2,1], aspect=1)
        viz!(ax3, m_floor;color=:lightgray)
        hidedecorations!(ax3)
        mcolor = vec(λ_infield_place)
        malpha = zeros(length(mcolor))
        malpha[isfinite.(mcolor)] .= 1.0
        mcolor[isnan.(mcolor)] .= 0.0
        viz!(ax3, m_floor, color=mcolor,alpha=malpha,colormap=colormap)
        Hippocampus.plot_pillars!(ax3)

        ax4 = Axis(lg3[3,1])
        Label(lg3[3,1, TopLeft()], "F")
        xx = [fill(1.0, length(X_place_nc));fill(2.0, length(pvc.λ_sub[spatial_field_idx,view_field_idx,:]))]
        yy = [X_place_nc;pvc.λ_sub[spatial_field_idx,view_field_idx,:]]
        cc = [fill(1, length(X_place_nc));fill(2, length(pvc.λ_sub[spatial_field_idx,view_field_idx,:]))]
        ss_pvc =  pvc.λ_covered[spatial_field_idx, view_field_idx] > percentile(pvc.λ_sub[spatial_field_idx, view_field_idx,:], 97.5)
        ss_pv =  μ_place_view > percentile(X_place_nc, 97.5)
        if ss_pvc
            marker_pvc = :star6
        else
            marker_pvc = :circle
        end
        if ss_pv
            marker_pv = :star6
        else
            marker_pv = :circle
        end
        if distr_plot == :boxplot
            boxplot!(ax4, xx,yy, show_outliers=false,orientation=:horizontal, color=cc, colormap=[:royalblue, :gray25],show_notch=true)
        else
            violin!(ax4, xx,yy, orientation=:horizontal, show_median=true, color=[:royalblue, :gray25][cc])
        end
        scatter!(ax4, [μ_place_view], [1.0],color=ccolors[spatial_field_idx], marker=marker_pv, label="Inside VF")
        scatter!(ax4, [μ_place], [1.0],color=:seagreen, label="Original VF")
        scatter!(ax4, [pvc.λ_covered[spatial_field_idx,view_field_idx]], [2.0], color=ccolors[view_field_idx], marker=marker_pvc,label="VF & PF")
        rowsize!(lg3, 3, 75)
        rowsize!(lg3,1, Relative(0.5))
        ax4.yticklabelsvisible = false
        ax4.yticksvisible = false
        ax4.xlabel = "Firing rate [Hz]"
        Hippocampus.link_cameras_lscene(fig)
        fig
    end
end

