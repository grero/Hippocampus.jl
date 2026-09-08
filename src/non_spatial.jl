using GLM
using StatsBase
using Random
using MultivariateStats


function get_poster_cue_response(celldirs::Vector{String};kwargs...)
    ncells = length(celldirs)
    nspikes = Vector{Vector{Int64}}(undef, ncells)
    posterid = Vector{Vector{Int64}}(undef, ncells)
    sessionid = fill(0, length(celldirs))
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    for (ii,celldir) in enumerate(celldirs)
        sessionid[ii] = findfirst(allsessiondirs[ii].==sessiondirs)
        _nspikes, _posterid, _outcome = get_poster_cue_response(celldir;kwargs...)
        nspikes[ii] = _nspikes[_outcome.==3]
        posterid[ii] = _posterid[_outcome.==3]
    end
    # if all sessiondirs are the same, just concatente
    if length(unique(sessionid)) == 1
        anspikes = cat(nspikes..., dims=2)
        aposterid = posterid[1]
    else
        anspikes = nspikes
        aposterid = posterid
    end
    anspikes, aposterid
end

function get_spikes(celldir::String;kwargs...)
    rp, udata = cd(DPHT.get_level_path("session", celldir)) do
        RippleData(), UnityData()
    end
    get_spikes(celldir, rp, udata;kwargs...)
end

function get_spikes(celldir::String, rp::RippleData, udata::UnityData;trial_start=2, trial_end=3)
    # get the spiketrains
    sp = cd(celldir) do
        sp = Spiketrain()
        sp
    end
    nt = numtrials(udata) 
    spiketimes = sp.timestamps/1000
    spikes = Vector{Vector{Float64}}(undef, nt)
    for i in 1:nt
        # TODO: Separate into correct and timeout trials here
        t0 = rp.timestamps[i,trial_start]
        t1 = rp.timestamps[i,trial_end]
        idx0 = searchsortedfirst(spiketimes, t0)
        idx1 = searchsortedlast(spiketimes, t1)
        spikes[i] = spiketimes[idx0:idx1] .- t0
    end
    spikes
end

function get_spikes(celldirs::Vector{String};kwargs...)
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    spikes = Vector{Vector{Float64}}[]
    for sessiondir in sessiondirs
        rp,udata = cd(sessiondir) do
            RippleData(), UnityData()
        end
        cidx = allsessiondirs.==sessiondir
        for celldir in celldirs[cidx]
            _spikes = get_spikes(celldir, rp, udata;kwargs...)
            push!(spikes, _spikes)
        end
    end
    spikes
end

function get_poster_cue_response(celldir::String;previous=false, future=0, tmin=0.0,kwargs...)
    # get the spiketrains
    sp, rp, udata = cd(celldir) do
        sp = Spiketrain()
        rp, udata = cd(DPHT.process_level("session")) do
            RippleData(), UnityData()
        end
        sp, rp, udata
    end
    nt = numtrials(udata)
    if previous
        i_start = 2
    else
        i_start = 1
    end
    if future > 0
        i_end = nt-future
    else
        i_end = nt
    end
    nspikes = zeros(Int64, i_end-i_start+1)
    posterid = zeros(Int64,i_end-i_start+1)
    outcome = zeros(Int64,i_end-i_start+1)
    spiketimes = sp.timestamps/1000
    # TODO: This doesn't work
    _posterid = rp.triggers[1,1] - 10
    for (j,i) in enumerate(i_start:i_end)
        # TODO: Separate into correct and timeout trials here
        t0,t1 = rp.timestamps[i,1:2]
        idx0 = searchsortedlast(spiketimes, t0)
        idx1 = searchsortedlast(spiketimes, t1)
        nspikes[j] = idx1-idx0
        if previous
            posterid[j] = _posterid 
            _posterid = rp.triggers[i,1]-10
        elseif future > 0
            posterid[j] = rp.triggers[i+future,1]-10
        else
            _posterid = rp.triggers[i,1]-10
            posterid[j] = _posterid 
        end
        outcome[j] = round(Int64,floor.(rp.triggers[i,3]/10))
    end
    nspikes, posterid, outcome
end

function compute_poster_selectivity(celldir::String;kwargs...)
    nspikes, posterid,outcome = get_poster_cue_response(celldir;kwargs...)
    cidx = outcome.==3
    vv0,vv = compute_poster_selectivity(nspikes[cidx],posterid[cidx];kwargs...)
end

function compute_poster_selectivity(nspikes::Vector{Int64}, posterid::Vector{Int64};nshuffles=100,kwargs...)
    y = indicatormat(posterid)
    model = glm(permutedims(y), nspikes, Poisson())
    vv0 = aic(model)
    vv = zeros(nshuffles)
    for i in 1:nshuffles
        _model = glm(permutedims(y), shuffle(nspikes), Poisson())
        vv[i] = aic(_model)
    end
    vv0, vv, model
end

function find_preferred_poster(nspikes::Vector{Int64}, posterid::Vector{Int64};k=1,nshuffles=1000)
    vv0,vv,model =  compute_poster_selectivity(nspikes, posterid;nshuffles=nshuffles)
    find_preferred_poster(model, nspikes, posterid;k=k,nshuffles=nshuffles)
end

function find_preferred_poster(model, nspikes::Vector{Int64}, posterid::Vector{Int64};k=1,nshuffles=1000)
    preferred = sortperm(model.pp.beta0,rev=true)[1:k]
    qlabel = fill(0, length(nspikes))
    tidx = findall(in(preferred[k]).(posterid))
    qlabel[tidx] .= 1
    # exlucde the 1:(k-1) preferred posters from further comparison
    aidx = findall((!in(preferred)).(posterid))
    qlabel[aidx] .= 2
    ww0,ww, model2 = compute_poster_selectivity(nspikes[qlabel.>0], qlabel[qlabel.>0];nshuffles=nshuffles)
    return ww0,ww, preferred, model2
end

function find_preferred_poster(celldir::String;kwargs...)
    pv_threshold = get(kwargs, :pv_threshold, 0.01)
    nspikes, posterid, outcome = get_poster_cue_response(celldir)
    cidx = outcome.==3
    vv0,vv,model = compute_poster_selectivity(nspikes[cidx], posterid[cidx];nshuffles=get(kwargs, :nshuffles, 1000))
    preferred = Int64[]
    if vv0 < percentile(vv, 100*pv_threshold)
        for k in 1:6
            ww0,ww,_preferred, m2 = find_preferred_poster(model, nspikes[cidx], posterid[cidx];k=k)
            if !(ww0 < percentile(ww, 100*pv_threshold))
                # return the previous result
                return vv0,vv,model, preferred
            else
                # update
                preferred = _preferred
            end
        end
    end
    return vv0, vv, model, Int64[]
end

function is_poster_selective(celldir::String;prctile=1, kwargs...)
    vv0,vv,model = compute_poster_selectivity(celldir;kwargs...)
    vv0 < percentile(vv, prctile), model
end

function is_poster_selective!(X::Vector{Bool}, celldirs::Vector{String};kwargs...)
    @showprogress "Finding selectivity..." for (i,celldir) in enumerate(celldirs)
        X[i] = is_poster_selective(celldir;kwargs...)
    end
    X
end

function plot_raster!(ax, spa::TrialAlignedSpiketrain, rp::RippleData;tmax=20.0, previous=false,future=false, markersize=5px, kwargs...)
    nt = numtrials(rp)
    posterid = rp.triggers[:,1] .- 10
    cue_period = spa.trigger_timestamps[:,2] .- spa.trigger_timestamps[:,1]
    trial_end = spa.trigger_timestamps[:,3] .- spa.trigger_timestamps[:,1]
    points = Tuple{Float64, Float64}[]
    color = to_colormap(:tab10)
    colors = typeof(color[1])[]
    if previous
        sidx = sortperm(collect(zip(posterid[1:end-1], trial_end[2:end])))
    elseif future
        sidx = sortperm(collect(zip(posterid[2:end], trial_end[1:end-1])))
    else
        sidx = sortperm(collect(zip(posterid, trial_end)))
    end
    vidx = invperm(sidx)
    #sposterid = posterid[sidx]
    if previous
        cue_period = cue_period[2:end][sidx]
    elseif future
        cue_period = cue_period[1:end-1][sidx]
    else
        cue_period = cue_period[sidx]
    end
    for (i,sp) in enumerate(spa.spiketimes)
        #only use correct trials
        if 3 .<= rp.triggers[i,3]/10 .< 4
            if trial_end[i] <= tmax
                for _sp in sp
                    if _sp <= tmax
                        if previous
                            if i > 1
                                push!(points, (_sp, vidx[i-1]))
                                push!(colors, color[posterid[i-1]])
                            end
                        elseif future
                            if i < nt
                                push!(points, (_sp, vidx[i]))
                                push!(colors, color[posterid[i+1]])
                            end
                        else
                            push!(points, (_sp, vidx[i]))
                            push!(colors, color[posterid[i]])
                        end
                    end
                end
            end
        end
    end
    if previous
        trial_end = trial_end[2:end][sidx]
    elseif future
        trial_end = trial_end[1:end-1][sidx]
    else
        trial_end = trial_end[sidx]
    end
    scatter!(ax, points,color=colors, markersize=get(kwargs, :markersize, 5px))
    scatter!(ax, Point2f.([(cp, i) for (i,cp) in enumerate(cue_period)]),color=:black, marker='|')
    #scatter!(ax, Point2f.([(trial_end[i], i) for i in findall(trial_end .<= tmax)]),color=:black, marker='|')
    if get(kwargs, :ylabelvisible, true)
        ax.ylabel = "Trialnr"
    end
    ax.yticklabelsvisible = false
end

function plot_psth!(ax, spa::TrialAlignedSpiketrain, rp::RippleData;tmax=20,binsize=0.05, window=5, previous=false, future=false,kwargs...)
    nt = numtrials(rp)
    weight,bins = compute_psth(spa, binsize, window)
    posterid = rp.triggers[:,1] .- 10
    color = to_colormap(:tab10)
    posterbins = unique(posterid)
    sort!(posterbins)
    idx1 = min(searchsortedlast(bins, tmax), size(weight,1))
    outcome = rp.triggers[:,3]
    cidx = findall(30 .<= outcome .< 40)
    for (i,b) in enumerate(posterbins)
        tidx = findall(posterid.==b)
        if previous
            # shift by one trial
            tidx .-= 1
            tidx = tidx[tidx .> 0]
        elseif future
            tidx .+= 1
            tidx = tidx[tidx .<= nt]
        end
        tidx = intersect(tidx, cidx) 
        μ = dropdims(mean(weight[:,tidx],dims=2),dims=2)
        μ ./= binsize*window
        lines!(ax, bins[1:idx1], μ[1:idx1];color=color[i])
    end
    if get(kwargs, :ylabelvisible, true)
        ax.ylabel = "Firing rate [Hz]" 
    end
end

function plot_poster_tuning!(ax, spa::TrialAlignedSpiketrain, rp::RippleData;tmin=0, tmax=:cue_onset, previous=false, future=false,kwargs...)
    nt = numtrials(rp)
    if (isa(tmax, Symbol) && tmax !== :cue_onset) || (isa(tmin, Symbol) && tmin !== :trial_end)
        if isa(tmax)
            error("Unknown symbol $tmax")
        else
            error("Unkonwn symbol $tmin")
        end
    end
    # just grab the spike counts during the cue period
    cue_onset = rp.timestamps[:,2] - rp.timestamps[:,1]
    trial_end = rp.timestamps[:,3] - rp.timestamps[:,1]
    spike_count = zeros(nt)
    for i in 1:nt
        if tmax == :cue_onset
            spike_count[i] = sum(tmin .< spa.spiketimes[i] .<= cue_onset[i])
        elseif tmin == :trial_end
            spike_count[i] = sum(trial_end[i] .< spa.spiketimes[i] .<= trial_end[i]+tmax)
        else
            spike_count[i] = sum(tmin .< spa.spiketimes[i] .<= tmax)
        end
    end
    posterid = rp.triggers[:,1] .- 10
    color = to_colormap(:tab10)
    posterbins = unique(posterid)
    sort!(posterbins)
    outcome = rp.triggers[:,3]
    cidx = findall(30 .<= outcome .< 40)
    label = fill(0, nt)
    for (i,b) in enumerate(posterbins)
        tidx = findall(posterid.==b)
        if previous
            # shift by one trial
            tidx .-= 1
            tidx = tidx[tidx .> 0]
        elseif future
            tidx .+= 1
            tidx = tidx[tidx .<= nt]
        end
        #tidx = intersect(tidx, cidx) 
        label[tidx] .= i
    end
    label[setdiff(1:nt, cidx)] .= 0
    # set up boxplot
    xx = label[label.>0]
    yy = spike_count[label.>0]
    boxplot!(ax, xx, yy, color=color[xx],show_notch=get(kwargs, :show_notch, false),
                        show_outliers=get(kwargs, :show_outliers, true))
end

function plot_raster_and_psth(args...;kwargs...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_raster_and_psth!(lg, args...;kwargs...)
        fig
    end
end

function plot_raster_and_psth!(lg::GridLayout, spa::TrialAlignedSpiketrain, rp::RippleData;show_psth=true, kwargs...)
    ax1 = Axis(lg[1,1])
    plot_raster!(ax1, spa, rp;kwargs...)
    ax2 = Axis(lg[2,1])
    if show_psth
        plot_psth!(ax2, spa, rp;kwargs...)
        linkxaxes!(ax1, ax2)
        if get(kwargs, :xlabelvisible, true)
            ax2.xlabel = "Time from cue [s]"
        end
    else
        plot_poster_tuning!(ax2, spa,rp;kwargs...)
    end
    ax1.xticklabelsvisible = false
    ax1.xticksvisible = false
    ax2.xticklabelsvisible = get(kwargs, :xticklabelsvisible, true)
end

function plot_raster_and_psth!(lg::GridLayout, spa1::TrialAlignedSpiketrain, spa2::TrialAlignedSpiketrain, rp::RippleData;show_psth=true, tmin=-1.0, tmax=20.0,show_connecting_lines=false,  kwargs...)
    ax11 = Axis(lg[1,1])
    ax12 = Axis(lg[1,2])
    # TODO: Plot both cue aligned and trial-end aligned raster

    plot_raster!(ax11, spa1, rp;tmin=tmin, tmax=tmax, kwargs...)
    plot_raster!(ax12, spa2, rp;tmin=-tmax, tmax=-tmin, kwargs...)
    if get(kwargs, :xlabelvisible, true)
        ax11.xlabel = "Time from cue [s]"
        ax12.xlabel = "Time from end [s]"
    end
    ax21 = Axis(lg[2,1],yticks=WilkinsonTicks(3))
    ax22 = Axis(lg[2,2], yticks=WilkinsonTicks(3))
    if show_psth
        plot_psth!(ax21, spa1, rp;kwargs...)
        plot_psth!(ax22, spa2, rp;kwargs...)
        linkxaxes!(ax11, ax21)
        linkxaxes!(ax12, ax22)
    else
        plot_poster_tuning!(ax21, spa1,rp;tmin=0.0, tmax=:cue_onset,kwargs...)
        plot_poster_tuning!(ax22, spa1,rp;tmin=:trial_end, tmax=1.0,kwargs...)
        if get(kwargs, :ylabelvisible, true)
            ax21.ylabel = "λ" 
        end
    end
    linkyaxes!(ax21,ax22)
    ax11.xticklabelsvisible = true 
    ax11.xticksvisible = true 
    ax12.xticklabelsvisible = true 
    ax12.xticksvisible = true 
    
    ax12.ylabelvisible = false
    ax22.ylabelvisible = false
    ax22.yticklabelsvisible = false
    ax21.xticklabelsvisible = get(kwargs, :xticklabelsvisible, true)
    ax22.xticklabelsvisible = get(kwargs, :xticklabelsvisible, true)
    ax21.xlabel = "Poster #"
    ax22.xlabel = "Poster #"

    #find the figure
    if isa(lg, Figure)
        parent_p = lg
    else
        parent_p = lg.parent
        while !isa(parent_p, Figure)
            parent_p = parent_p.parent
        end
    end
    if show_connecting_lines
        ll = ax21.finallimits[]
        pos_fig_1 = pos_fig_obs(ax11, 0.0, 0.0)
        pos_fig_2 = pos_fig_obs(ax21, ll.origin[1], ll.origin[2]+ll.widths[2])
        pts1 = @lift [$pos_fig_1, $pos_fig_2]
        lw = 1.0
        #lines!(lg.parent.parent.parent.scene, pts1,color=:black, linewidth=lw)
        lines!(parent_p.scene, pts1,color=:black, linewidth=lw)

        pos_fig_3 = pos_fig_obs(ax11, 1.0, 0.0)
        pos_fig_4 = pos_fig_obs(ax21, ll.origin[1]+ll.widths[1], ll.origin[2]+ll.widths[2])
        pts2 = @lift [$pos_fig_3, $pos_fig_4]
        #lines!(lg.parent.parent.parent.scene, pts2,color=:black, linewidth=lw)
        lines!(parent_p.scene, pts2,color=:black, linewidth=lw)

        ll = ax22.finallimits[]
        pos_fig_5 = pos_fig_obs(ax12, 0.0, 0.0)
        pos_fig_6 = pos_fig_obs(ax22, ll.origin[1], ll.origin[2]+ll.widths[2])
        pts3 = @lift [$pos_fig_5, $pos_fig_6]
        #lines!(lg.parent.parent.parent.scene, pts3,color=:black, linewidth=lw)
        lines!(parent_p.scene, pts3,color=:black, linewidth=lw)

        pts4 = @lift begin
            ll = $(ax22.finallimits)
            # I'm not sure why this doens't work some time
            # I'm just going to hard code 
            pos_fig_7 = $(pos_fig_obs(ax12, 1.0, 0.0))
            # This does not work
            #pos_fig_8 = $(pos_fig_obs(ax22, ll.origin[1]+ll.widths[1], ll.origin[2]+ll.widths[2]))
            # for some reason, some times ll.origin[1] + ll.widths[1] == 10.0
            pos_fig_8 = $(pos_fig_obs(ax22, 6.9, ll.origin[2]+ll.widths[2]))
            #@show ll.origin[1] + ll.widths[1]
            [Point2f(pos_fig_7), Point2f(pos_fig_8)]
        end
        #lines!(lg.parent.parent.parent.scene, pts4,color=:black, linewidth=lw)
        lines!(parent_p.scene, pts4,color=:black, linewidth=lw)
    else
        # just show a gray bar
        linesegments!(ax11, [(Point2f(0.0, -6.0), Point2f(1.0, -6.0))], color=:gray, linewidth=5.0)
    end
end

function plot_raster_and_psth!(lg::GridLayout, celldir::String;kwargs...)
    spa1,spa2, rp = cd(celldir) do
        rp = cd(DPHT.process_level("session")) do
            RippleData()
        end
        spa1 = TrialAlignedSpiketrain(;alignto=1,Δt1=1.0)
        spa2 = TrialAlignedSpiketrain(;alignto=3)
    spa1, spa2, rp
    end
    plot_raster_and_psth!(lg, spa1, spa2, rp;kwargs...)
end

function population_decoder(nspikes::Matrix{T}, posterid::Vector{Int64};nruns=1) where T <: Real
    nt = length(posterid)
    nc = maximum(posterid)
    # whiten
    pca = fit(PCA, nspikes)
    y = predict(pca, nspikes)
    ncells = size(nspikes,1)

    ntrain = round(Int64, 0.8*nt)
    ntest = nt - ntrain
    perf = fill(0.0, nruns)
    w = fill(0.0, ncells, nruns)
    for r in 1:nruns
        train_idx = shuffle(1:nt)[1:ntrain]
        sort!(train_idx)
        test_idx = setdiff(1:nt, train_idx)
        lda = fit(MulticlassLDA, nc, y[:,train_idx], posterid[train_idx])
        _w = pca.proj*lda.proj
        w[:,r] .= dropdims(sum(abs2,_w,dims=2),dims=2)
        w[:,r] ./= sum(w[:,r],dims=1)
        z = predict(lda, y[:,test_idx])
        cmeans = predict(lda, classmeans(lda))
        for j in 1:ntest
            d = dropdims(sum(abs2, z[:,j] .- cmeans, dims=1),dims=1)
            k = argmin(d)
            perf[r] += k == posterid[test_idx[j]]
        end
        perf[r] /= ntest
    end
    perf,w
end

function population_decoder(allcelldirs::Vector{String};nruns=10, kwargs...)
    # find the unique sessions
    allsessiondirs = DPHT.get_level_path.("session", allcelldirs)
    sessiondirs = unique(allsessiondirs)
    nsessions = length(sessiondirs)
    # decode each session separately
    ncells = fill(0, nsessions) 
    perf = fill(0.0, nruns, nsessions)
    nt = fill(0, nsessions)
    coding_contrib = fill(0.0, length(allcelldirs),nruns)
    for (ii,sessiondir) in enumerate(sessiondirs)
        cidx = findall(allsessiondirs.==sessiondir)
        ncells[ii] = length(cidx)
        nspikes, posterid = get_poster_cue_response(allcelldirs[cidx];kwargs...)
        nt[ii] = size(nspikes,1)
        perf[:,ii],w = population_decoder(permutedims(nspikes), posterid;nruns=nruns)
        # compute relative contribution for each cell
        coding_contrib[cidx,:] .= w
    end
    perf, nt, ncells, coding_contrib, sessiondirs
end

function population_decoder(nspikes::Vector{Vector{Int64}}, posterid::Vector{Vector{Int64}})

end