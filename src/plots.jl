using CairoMakie


function plot_place_field(celldir::String;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(1000,250))
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
    ax0 = Axis(lg[1,1])
    boxplot!(ax0, fill(1.0, length(sic.sic)), sic.sic;show_outliers=false, show_notch=true,color=:gray)
    ax0.xticklabelsvisible = false
    ax0.xticksvisible = false
    ax0.bottomspinevisible = false
    scatter!(ax0, [1.0],[sic.sic0],color=:red)
    ax0.ylabel = "SIC"

    ax1 = Axis(lg[1,2])
    ax2 = Axis(lg[1,4])
    ax3 = Axis(lg[1,6])

    Z = sm.weight./sm.occupancy
    viz!(ax1, mm;showsegments=false, color=:lightgray)
    viz!(ax1, mm;showsegments=false, color=Z)
    Colorbar(lg[1,3], colorrange=(extrema(filter(isfinite, Z))),label="Firing rate [Hz]")
    ax1.title = "Raw"
    Zg = smg.weight./smg.occupancy
    Zg[smg.unvisited] .= NaN 
    viz!(ax2, mm;showsegments=false, color=:lightgray)
    viz!(ax2, mm;showsegments=false, color=Zg)
    Colorbar(lg[1,5], colorrange=(extrema(filter(isfinite, Zg))), label="Firing rate [Hz]")

    ax2.title = "σ = $(σ)"
    Za = sma.weight./sma.occupancy
    Za[sma.unvisited] .= NaN
    viz!(ax3, mm;showsegments=false, color=:lightgray)
    viz!(ax3, mm;showsegments=false, color=Za)
    Colorbar(lg[1,7], colorrange=(extrema(filter(isfinite, Za))), label="Firing rate [Hz]")
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
    colsize!(lg, 1, Relative(0.1))
    ax1,ax2,ax3
end

"""
Plot a summary of the spatial selectivity across all cells

One panel showing the distribution of SIC scores
Some example cells. One near the top, i.e. the most spatially selective, one near the median and one near the bottom?
"""
function plot_spatial_summary!(lg,celldirs::Union{Vector{String},Nothing}=nothing;kwargs...)
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
    if isfile(fname)
        res,sic = JLD2.load(fname, "res","sic")
    else
        res = issignificant(SpatialInformationContent, celldirs;skip_error=true, kwargs...)
        sic = get_sic(SpatialInformationContent, celldirs;skip_error=true, kwargs...)
        JLD2.save(fname, Dict("res"=>res, "sic"=>sic))
    end

    lg1 = GridLayout(lg[1,1])
    Label(lg1[1,1,TopLeft()], "A")
    ax1 = Axis(lg1[1,1])
    _colors = parse.(Colorant, [:grey, :red])
    cc = fill(_colors[1], length(sic))
    cc[res] .= _colors[2] 
    xx = rand(length(cc))
    yy = sic
    # TODO: Plot this differently
    scatter!(ax1, xx, yy, color=cc)
    Legend(lg1[1,2], [MarkerElement(marker=:circle, color=q) for q in _colors],["Non-sig","sig"],tellwidth=false, tellheight=false,
           valign=:top, halign=:right, framevisible=true, padding=(10.0, 10.0, 10.0, 10.0))
    #rainclouds!(ax1, xx, yy)
    #ax1.xticks = (1:2, ["Non-sig", "sig"])
    ax1.ylabel = "SIC"
    ax1.xticklabelsvisible = false
    ax1.xticksvisible = false
    ax1.bottomspinevisible = false
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
    plot_place_field!(lg2, celldirs[cidx1];ylabel="Cell 1")

    lg3 = GridLayout(lg[3,1])
    cidx2 = non_place_idx[argmax(sic[non_place_idx])]
    plot_place_field!(lg3, celldirs[cidx2];ylabel="Cell 2")
    lg4 = GridLayout(lg[4,1])
    # median
    cidx3 = place_idx[sortperm(sic[place_idx])[div(length(place_idx),2)]]
    plot_place_field!(lg4, celldirs[cidx3];ylabel="Cell 3")

    # TODO: Indicate these two cells in the scatter plot
    cidx = [cidx1, cidx2, cidx3]
    scatter!(ax1, xx[cidx], yy[cidx], markersize=10px, strokecolor=:black, color=cc[cidx], strokewidth=2.0)
    annotation!(ax1, xx[cidx], yy[cidx]; text=["1","2","3"]) 
end

function plot_spatial_summary(celldirs::Union{Vector{String},Nothing}=nothing;kwargs...)
    with_theme(plot_theme) do
        fig = Figure(size=(1035, 1034))
        lg = GridLayout(fig[1,1])
        plot_spatial_summary!(lg, celldirs;kwargs...)
        fig
    end
end
