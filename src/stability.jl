using OptimalTransport

struct SpatialMapStability
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    ccm::Float64
    ccms::Vector{Float64}
end

function get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:all, smooth=false, α=0.5, niter=50,kwargs...)
    jm1 = JointMap(vpvrp, jocc,jocc_filtered;use_trials=use_trials)
    spm1 = Hippocampus.SpatialMapNew(jm1,m_floor)
    if smooth
        spml = SmoothedMap(spm1;method=:laplace, α=α, niter=niter)
        λ1 = get_rate_map(spml)
    else
        λ1 = get_rate_map(spm1)
    end
    λ1
end

function get_cross_correlation(λ1, λ2)
    qidx = (isfinite.(λ1)).&(isfinite.(λ2))
    cc = crosscor(λ1[qidx], λ2[qidx])
    cm = maximum(cc)
end

function compare_maps(λ1, λ2;normalize=true)
    fidx1 = isfinite.(λ1)
    fidx2 = isfinite.(λ2)
    fidx = fidx1.&fidx2
    if normalize
        λ1n =  (λ1.-mean(λ1[fidx1]))
        λ1n ./=norm(λ1n[fidx1])
        λ2n = (λ2.-mean(λ2[fidx2]))
        λ2n ./= norm(λ2n[fidx2])
    else
        λ1n = λ1
        λ2n = λ2
    end
    λ1n[fidx]'*λ2n[fidx]
end

function get_cross_correlation_1(λ1::AbstractVector{<:Real}, λ2::AbstractVector{<:Real}, mm::SimpleMesh;max_x_shift=nothing, max_y_shift=nothing)
    nn = nelements(mm)
    max_shift = div(round(Int64, ceil(sqrt(nn))),2)
    if max_x_shift === nothing
        max_x_shift = max_shift
    end
    if max_y_shift === nothing
        max_y_shift = max_shift
    end
    # do cross-correlation the brute force way, by explicitly shifting one map with respect to the other
    nn = nelements(mm)

    v = get_tangential_space(mm);
    A = adjacencymatrix(mm)
    λ1s = copy(λ1)
    qq = zeros(2*max_y_shift+1, 2*max_x_shift+1)
    shifts = Matrix{Tuple{Float64, Float64}}(undef, 2*max_y_shift+1, 2*max_x_shift+1)
    qq[max_y_shift+1, max_x_shift+1] = compare_maps(λ1, λ2) 
    shifts[max_y_shift+1, max_x_shift+1] = (0.0, 0.0)
    for (j,sx) in enumerate(1:1:max_x_shift)
        λ1s = shift_map(λ1, Vec2(-1.0, 0.0), v, A;nsteps=sx)
        qq[max_y_shift+1,j] = compare_maps(λ1s, λ2)
        shifts[max_y_shift+1,j] = (0.0, -sx)
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, -1.0), v, A;nsteps=sy)
            qq[i,j] = compare_maps(λ1ss, λ2)
            shifts[i,j] = (-sy,-sx)
        end
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, 1.0), v, A;nsteps=sy)
            qq[max_y_shift+i+1,j] = compare_maps(λ1ss, λ2)
            shifts[max_y_shift+i+1,j] = (-sx, sy)
        end
    end
     for (j,sx) in enumerate(1:1:max_x_shift)
        λ1s = shift_map(λ1, Vec2(1.0, 0.0), v, A;nsteps=sx)
        qq[max_y_shift+1,max_x_shift+j+1] = compare_maps(λ1s, λ2)
        shifts[max_y_shift+1,max_x_shift+j+1] = (0.0, sx)
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, -1.0), v, A;nsteps=sy)
            qq[i,max_x_shift+j+1] = compare_maps(λ1ss, λ2)
            shifts[i,max_x_shift+j+1] = (-sy,sx)
        end
        for (i,sy) in enumerate(1:1:max_y_shift)
            λ1ss = shift_map(λ1s, Vec2(0.0, 1.0), v, A;nsteps=sy)
            qq[max_y_shift+i+1,max_x_shift+j+1] = compare_maps(λ1ss, λ2)
            shifts[max_y_shift+i+1,max_x_shift+j+1] = (sy,sx)
        end
    end
    qq,shifts
end

function get_cross_correlation(λ1::AbstractVector{<:Real}, λ2::AbstractVector{<:Real}, mm::SimpleMesh)
    # get the bins
    cp = Point2f.(Tuple.(centroid.(mm)))
    # get the binsize
    binsize = [v.val for v in sqrt.(measure.(mm))]

    D = maximum(norm.(cp .- permutedims(cp)))
    max_shift = round(Int64,ceil(floor(D/sqrt(2)/2)/binsize[1]))
    qq = zeros(2*max_shift+1, 2*max_shift+1)
    lags = Matrix{Tuple{Float64, Float64}}(undef, size(qq)...)
    _binsize = binsize[1]
    for (i,sx) in enumerate(-max_shift:max_shift)
        λ1s,_ = shift_mass(λ1, cp, binsize, Vec2(sx*_binsize, 0.0)) 
        for (j,sy) in enumerate(-max_shift:max_shift)
            λ1ss,_ = shift_mass(λ1s, cp, binsize, Vec2(0.0, sy*_binsize))
            qq[j,i] = compare_maps(λ1ss, λ2) 
            lags[j,i] = (sy*_binsize,sx*_binsize)
        end
    end
    qq,lags
end

"""
    get_dissimilarity(λ1, λ2, mm::SimpleMesh)

Compute the dissimilarity between maps `λ1` and `λ2` defined on the mesh `mm`, using the sinkhole algorithm.
"""
function get_dissimilarity(λ1, λ2, mm::SimpleMesh)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    D = Hippocampus.distancematrix(mm;between_centroids=false)
    get_dissimilarity(λ1, λ2,D[fidx,fidx])
end

function get_dissimilarity(λ1, λ2, D::AbstractMatrix{<:Real})
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    f = λ1[fidx]./norm(λ1[fidx])
    g = λ2[fidx]./norm(λ2[fidx])
    P = sinkhorn_unbalanced(f,g,D, 1.0, 1.0, 0.01)
    #P = sinkhorn(f,g,D, 0.01)
    P,D
end

function get_null_distr(λ1, λ2, mm::SimpleMesh;nruns=1000)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    D = Hippocampus.distancematrix(mm;between_centroids=false,idx=fidx)
    λ1s = copy(λ1)
    λ2s = copy(λ2)
    fidx = findall((isfinite.(λ1)).&(isfinite.(λ2)))
    dd = zeros(nruns)
    for r in 1:nruns
        for i in fidx
            if rand() < 0.5
                λ1s[i] = λ2[i]
            else
                λ1s[i] = λ1[i]
            end
            if rand() < 0.5
                λ2s[i] = λ1[i]
            else
                λ2s[i] = λ2[i]
            end
        end
        P,D = get_dissimilarity(λ1s, λ2s, D)
        dd[r] = sum(P.*D)
    end
    dd
end

function process_kwargs(::Type{SpatialMapStability},h::UInt32=zero(UInt32); nshuffles=1000, smooth=false, α=0.1, niter=50, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles => nshuffles),h)
    h = crc32c(string(:smooth=>smooth),h)
    if smooth
        h = crc32c(string(:α=>α),h)
        h = crc32c(string(:niter=>niter),h)
    end
    h
end

function SpatialMapStability(;redo=fname->false, do_save=true, nshuffles=1000, kwargs...)
    fname = "spatial_map_stability.jld2"
    h = process_kwargs(SpatialMapStability;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStability, fname)
    else
        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
        jocc,qdata = cd(DPHT.process_level("session")) do
            jocc = JointOccupancy(;kwargs...)
            qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
            jocc,qdata
        end
        jocc_filtered = JointFilteredOccupancy(jocc, qdata;kwargs...)

        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        obj = SpatialMapStability(vpvrp, jocc, jocc_filtered, m_floor;nshuffles=nshuffles, kwargs...)
        if do_save
            save_jld2(obj, fname;nshuffles=nshuffles, kwargs...)
        end
    end
    obj
end

function get_spatial_map_stability(;nshuffles=10_000, kwargs...)
    # FIXME: There is something fishy with this function that prevents me from getting the result in the REPL
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    trial_start= get(kwargs, :trial_start, 2)
    # assume we are in the correct directory; load whatever we need
    sp = Spiketrain()

    rp = cd(DPHT.process_level(level(RippleData))) do
        RippleData()
    end
    rs1 = RandomlyShiftedSpiketrains(;nshifts=nshuffles, use_trials = :firstHalf, trial_start=trial_start)
    rs2 = RandomlyShiftedSpiketrains(;nshifts=nshuffles, use_trials = :secondHalf, trial_start=trial_start)

    unity_gaze_data = cd(DPHT.process_level("session")) do
        UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
    end
    jocc = cd(DPHT.process_level(level(JointOccupancy))) do
        JointOccupancy(;redo=fname->false, kwargs...)
    end
    jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
    vpvrp = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...) 

    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    P,D = get_dissimilarity(λ1, λ2, m_floor)
    ccm = sum(P.*D)
    ccms = zeros(nshuffles)
    for (ii,(sp1,sp2)) in enumerate(zip(eachcol(rs1.timestamps), eachcol(rs2.timestamps)))
        vpvrp1 = ViewAndPlaceRepresentationNew(sp1/1000.0,rp,unity_gaze_data;kwargs...) 
        λ1s = get_spatial_map(vpvrp1, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
        vpvrp2 = ViewAndPlaceRepresentationNew(sp2/1000.0,rp,unity_gaze_data;kwargs...) 
        λ2s = get_spatial_map(vpvrp2, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
        P,D = get_dissimilarity(λ1s, λ2s, m_floor)
        ccms[ii] = sum(P.*D)
    end
    ccm, ccms
end

function get_spatial_map_stability2(;nshuffles=10_000, kwargs...)
    # FIXME: There is something fishy with this function that prevents me from getting the result in the REPL
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    jm1 =  JointMap(;use_trials=:firstHalf, kwargs...)
    spm1 = SpatialMapNew(jm1,m_floor);
    λ1 = get_rate_map(spm1);
    jm2 =  JointMap(;use_trials=:secondHalf, kwargs...)
    spm2 = SpatialMapNew(jm2,m_floor); 
    λ2 = get_rate_map(spm2);

    P,D = get_dissimilarity(λ1, λ2, m_floor)
    ccm = sum(P.*D)
    qq, lags = get_cross_correlation(λ1, λ2, m_floor)
    ccq = maximum(qq)
    ccqs = zeros(nshuffles)
    ccms = zeros(nshuffles)
    for ii in 1:nshuffles
        #λ1s = shuffle(λ1)
        λ2s = shuffle(λ2)
        P,D = get_dissimilarity(λ1, λ2s, m_floor)
        ccms[ii] = sum(P.*D)
        qq, _ = get_cross_correlation(λ1, λ2s, m_floor)
        ccqs[ii] = maximum(qq)
    end
    ccm, ccms, ccq, ccqs
end


struct SpatialMapStabilityNew
    λ1::Vector{Float64}
    λ2::Vector{Float64}
    cc_ff::Vector{Float64}
    cc_fs::Vector{Float64}
    ds_ff::Vector{Float64}
    ds_fs::Vector{Float64}
end

function process_kwargs(::Type{SpatialMapStabilityNew}, h::UInt32=zero(UInt32);nshuffles=10_000, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    h
end


function SpatialMapStabilityNew(;redo=fname->false, do_save=true, nshuffles=10_000, kwargs...)
    fname = "spatial_map_stability_new.jld2"
    h = process_kwargs(SpatialMapStabilityNew;nshuffles=nshuffles, kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(SpatialMapStabilityNew, fname)
    else

        nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
        m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))

        sp = Spiketrain()

        rp = cd(DPHT.process_level(level(RippleData))) do
            RippleData()
        end

        unity_gaze_data = cd(DPHT.process_level("session")) do
            UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
        end
        jocc = cd(DPHT.process_level(level(JointOccupancy))) do
            JointOccupancy(;redo=fname->false, kwargs...)
        end
        nt = length(jocc.index)
        nthalf = div(nt,2)
        trialidx = [1:nt;]
        trialidx = collect(1:nt)
        jocc_filtered = JointFilteredOccupancy(jocc, unity_gaze_data;kwargs...)
        vpvrp = ViewAndPlaceRepresentationNew(sp,rp,unity_gaze_data;kwargs...) 

        cc_ff = zeros(nshuffles)
        ds_ff = zeros(nshuffles)
        cc_fs = zeros(nshuffles)
        ds_fs = zeros(nshuffles)
        @showprogress for i in 1:nshuffles
            # random maps from the first half of the trials
            trialidx1 = sort(shuffle(1:nthalf)[1:div(nthalf,2)])
            trialidx2 = setdiff(1:nthalf, trialidx1)
            λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx1,kwargs...)
            λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)
            P,D = get_dissimilarity(λ1, λ2, m_floor)
            ds_ff[i] = sum(P.*D)
            qq, lags = get_cross_correlation(λ1, λ2, m_floor)
            cc_ff[i] = maximum(qq)

            #second half
            trialidx2 = sort(shuffle(nthalf:1+nt)[1:div(nt-nthalf+1,2)])
            λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)

            P,D = get_dissimilarity(λ1, λ2, m_floor)
            ds_fs[i] = sum(P.*D)
            qq, lags = get_cross_correlation(λ1, λ2, m_floor)
            cc_fs[i] = maximum(qq)
        end
        λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
        λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
        obj = SpatialMapStabilityNew(λ1, λ2,cc_ff, cc_fs, ds_ff, ds_fs)
        if do_save
            save_jld2(obj, fname;kwargs...)
        end
    end
    obj
end

function get_spatial_map_stability(vpvrp, jocc, jocc_filtered, m_floor::SimpleMesh;kwargs...)
    # compute stability of mean spatial response for 1st and 2nd half of the trials, compare with equivalent number of random trials
    # If there is a systematic shift 
    nt = length(jocc.index)
    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    P,D = get_dissimilarity(λ1, λ2, m_floor)
    ccm = sum(P.D)
end

function SpatialMapStability(vpvrp, jocc, jocc_filtered, m_floor::SimpleMesh;nshuffles=1000, kwargs...)
    # compute stability of mean spatial response for 1st and 2nd half of the trials, compare with equivalent number of random trials
    # If there is a systematic shift 
    nrefinements = get(kwargs, :nrefinements, (p=3,g=2))
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    nt = length(jocc.index)
    λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:firstHalf,kwargs...)
    λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=:secondHalf,kwargs...)
    qq,lags = get_cross_correlation(λ1, λ2, m_floor)
    ccm = maximum(qq)
    ccms = zeros(nshuffles)
    nthalf = div(nt,2)
    trialidx = [1:nt;]
    @showprogress for i in 1:length(ccms)
        #trialidx1 = sort(shuffle(trialidx)[1:nthalf])
        #trialidx2 = setdiff(trialidx, trialidx1)
        #λ1 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx1,kwargs...)
        #λ2 = get_spatial_map(vpvrp, jocc, jocc_filtered,m_floor;use_trials=trialidx2,kwargs...)
        qqs,lags = get_cross_correlation(shuffle(λ1), shuffle(λ2),m_floor)
        ccms[i] = maximum(qqs)
    end
    SpatialMapStability(λ1, λ2, ccm, ccms)
end

## plots

function plot_stability_summary(celldirs::Vector{String};_plot_theme=plot_theme, kwargs...)
    kargs = (nrefinements=(p=1,g=1), min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, smooth=false) 
    spm_stability = map(celldirs) do celldir
        spm_stability = cd(celldir) do
            Hippocampus.SpatialMapStability(;kargs...)
        end
        spm_stability.ccm
    end
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=1))
    with_theme(_plot_theme) do
        fig = Figure(size=(1.5*650,1.5*300))
        ax = Axis(fig[1,1])
        hist!(ax, spm_stability;color=:gray)
        ax.xlabel = "Stability"
        # show example of the 5th percentile, the median, and the 95th percentile cells in terms of stability
        lg = GridLayout(fig[1,2])
        lg1 = GridLayout(lg[1,1], alignmode=Outside(5))
        lg2 = GridLayout(lg[1,2], alignmode=Outside(5))
        lg3 = GridLayout(lg[1,3], alignmode=Outside(5))
        idx0 = argmin(norm.(spm_stability .- percentile(spm_stability, 5)))
        idx1 = argmin(norm.(spm_stability .- percentile(spm_stability, 50)))
        idx2 = argmin(norm.(spm_stability .- percentile(spm_stability, 95)))
        @show idx0 idx1 idx2
        # indicate these points on the histogram
        colors = [:pink, :red, :orange]
        vlines!(ax, spm_stability[[idx0,idx1,idx2]], color=colors)
        # draw boxes around the corresponding plots
        Makie.Box(lg[1,1], color=(:white, 0.0), strokecolor=colors[1])
        Makie.Box(lg[1,2], color=(:white,0.0), strokecolor=colors[2])
        Makie.Box(lg[1,3], color=(:white,0.0), strokecolor=colors[3])
        for (ii,(idx,_lg)) in enumerate(zip([idx0, idx1, idx2],[lg1,lg2,lg3]))
            jm1,jm2 = cd(celldirs[idx])  do
                jm1 = JointMap(;use_trials=:firstHalf, kargs...)
                jm2 = JointMap(;use_trials=:secondHalf, kargs...)
                jm1, jm2
            end
            spm1 = SpatialMapNew(jm1, m_floor)
            spm2 = SpatialMapNew(jm2, m_floor)
            ax1 = Axis(_lg[1,1], aspect=1.0)
            ax2 = Axis(_lg[2,1], aspect=1.0)
            for _ax in [ax1, ax2]
                hidedecorations!(_ax)
                _ax.bottomspinevisible = false
                _ax.leftspinevisible = false
            end
            λ1 = get_rate_map(spm1)
            λ2 = get_rate_map(spm2)
            cr = extrema(filter(isfinite, [λ1;λ2]))
            viz!(ax1, m_floor;color=:darkgray)
            viz!(ax1, m_floor;color=get_rate_map(spm1),colormap=:rain, colorrange=cr)
            viz!(ax2, m_floor;color=:darkgray)
            viz!(ax2, m_floor;color=get_rate_map(spm2),colormap=:rain, colorrange=cr)
        end
        lgf = GridLayout(fig[1,3])
        Label(lgf[1,1], "First half", rotation=π/2, tellheight=false)
        Label(lgf[2,1], "Second half", rotation=π/2, tellheight=false)
        colsize!(fig.layout, 1, Relative(0.4))
        fig
    end
end


function plot_stability(args...;_plot_theme=plot_theme,kwargs...)
    with_theme(_plot_theme) do
        fig = Figure(size=(600,400))
        lg = GridLayout(fig[1,1])
        plot_stability!(lg, args...;kwargs...)
        fig
    end
end

function plot_stability!(lg, vpvrp::ViewAndPlaceRepresentationNew, jocc, jocc_filtered, m_floor::SimpleMesh;kwargs...)
    λ1 = Hippocampus.get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:firstHalf)
    λ2 = Hippocampus.get_spatial_map(vpvrp, jocc, jocc_filtered, m_floor;use_trials=:secondHalf)
    plot_stability!(lg, λ1, λ2, m_floor;kwargs...)
end

function plot_stability!(lg, spm::SpatialMapStabilityNew, m_floor::SimpleMesh;plot_cross_correlation=true)
    if plot_cross_correlation
        (s1,s2) = (spm.cc_ff, spm.cc_fs)
    else
        (s1,s2) = (spm.ds_ff, spm.ds_fs)
    end
    plot_stability!(lg, spm.λ1, spm.λ2, m_floor, s1,s2)
end

function plot_stability!(lg, λ1::AbstractVector{<:Real},λ2::AbstractVector{<:Real}, m_floor, s1::AbstractVector{<:Real}, s2::AbstractVector{<:Real};ylabel="Cross-correlation", kwargs...)
    lg1 = GridLayout(lg[1,1])
    plot_stability!(lg1, λ1, λ2, m_floor;kwargs...)
    ax = Axis(lg1[2,1:2])
    Label(lg1[1,1,Top()], "First half")
    Label(lg1[1,2,Top()], "Second half")
    xx = [fill(1.0, length(s1));fill(2.0, length(s2))]
    yy = [s1;s2]
    fidx = isfinite.(yy)
    rainclouds!(ax, xx[fidx], yy[fidx])
    ax.ylabel = ylabel
    ax.xticklabelsvisible = true 
    ax.xticksvisible =  false 
    ax.bottomspinevisible = false 
    ax.xticks = ([1,2], ["fh→fh","fh→sh"])
    Label(lg1[1,1,TopLeft()], "A")
    Label(lg1[2,1,TopLeft()], "B")
end

function plot_stability!(lg, λ1::AbstractVector{<:Real},λ2::AbstractVector{<:Real}, m_floor;kwargs...)
    cr = extrema(filter(isfinite, [λ1;λ2]))
    colormap = get(kwargs, :colormap, :rain)
    ax1 = Axis(lg[1,1], aspect=1.0)
    viz!(ax1, m_floor;color=:darkgray)
    viz!(ax1, m_floor;color=λ1, colormap=colormap, colorrange=cr)

    ax2 = Axis(lg[1,2], aspect=1.0)
    viz!(ax2, m_floor;color=:darkgray)
    viz!(ax2, m_floor;color=λ2, colormap=colormap, colorrange=cr)

    Colorbar(lg[1,3], colorrange=cr, colormap=colormap, label="Firing rate [Hz]")
    for _ax in [ax1, ax2]
        hidedecorations!(_ax)
        _ax.bottomspinevisible = false
        _ax.leftspinevisible = false
    end
end

function plot_stability(celldir::String;kargs...)
    
end