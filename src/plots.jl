using CairoMakie


function plot_place_field(celldir::String;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_place_field!(lg, celldir;kwargs...)
        fig
    end
end

"""
Illustrate place fields for the cell in `celldir`.

One panel with the raw firing rate map, one with a gaussian smoothed firing rate map and one with the adaptively smoothed firing rate map
"""
function plot_place_field!(lg, celldir::String;nrefinements=(p=3,g=2),σ=3, α=1000.0.^2, ylabel="", kwargs...)
    sm, smg, sma,mm, sic = cd(celldir) do
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
        sic = compute_skaggs_sic(SpatialInformationContent,10_000;kwargs...)
        mm = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        #vpvrpb = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...)
        jmb = JointMap(;kwargs...)
        # TODO: I need to be able to load this directly
        spmb = SpatialMapNew(jmb,mm)
        smg = SmoothedMap(spmb;method=:gaussian, σ=σ)
        sma = SmoothedMap(spmb;method=:adaptive, α=α)
        spmb,smg,sma,mm, sic
    end
    @show sic.sic0 percentile(sic.sic, 95)
    # add SIC with distribution
    ax1 = Axis(lg[1,1])
    ax2 = Axis(lg[1,3])
    ax3 = Axis(lg[1,5])

    Z = sm.weight./sm.occupancy
    viz!(ax1, mm;showsegments=false, color=:lightgray)
    viz!(ax1, mm;showsegments=false, color=Z)
    Colorbar(lg[1,2], colorrange=(extrema(filter(isfinite, Z))),label="Firing rate [Hz]")
    ax1.title = "Raw"
    Zg = smg.weight./smg.occupancy
    Zg[smg.unvisited] .= NaN 
    viz!(ax2, mm;showsegments=false, color=:lightgray)
    viz!(ax2, mm;showsegments=false, color=Zg)
    Colorbar(lg[1,4], colorrange=(extrema(filter(isfinite, Zg))), label="Firing rate [Hz]")

    ax2.title = "σ = $(σ)"
    Za = sma.weight./sma.occupancy
    Za[sma.unvisited] .= NaN
    viz!(ax3, mm;showsegments=false, color=:lightgray)
    viz!(ax3, mm;showsegments=false, color=Za)
    Colorbar(lg[1,6], colorrange=(extrema(filter(isfinite, Za))), label="Firing rate [Hz]")
    ax3.title = "α = $α"
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
    ax1,ax2,ax3
end