using CrossTemporalDecoding
"""
    get_data(celldirs::Vector{String};kwargs...)

Get spiketrains and poster combinations for all celldirs
"""
function get_data(celldirs::Vector{String};correct_after_correct_only=true,kwargs...)
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    spikes = Vector{Vector{Vector{Vector{Float64}}}}(undef, length(sessiondirs))
    poster_labels = Vector{Vector{Tuple{Int64, Int64}}}(undef, length(sessiondirs))
    for (ii,sessiondir) in enumerate(sessiondirs)
        rp,udata = cd(sessiondir) do
            RippleData(), UnityData()
        end
        _poster_labels,correct_trial_idx = get_poster_combinations(udata)
        if correct_after_correct_only
            correct_trial_idx = findall(correct_trial_idx[2:end].&correct_trial_idx[1:end-1]) .+ 1
        end
        poster_labels[ii] = _poster_labels[correct_trial_idx]
        cidx = allsessiondirs.==sessiondir
        spikes[ii] = Vector{Vector{Float64}}(undef, sum(cidx))
        for (jj,celldir) in enumerate(celldirs[cidx])
            _spikes = get_spikes(celldir;kwargs...)
            spikes[ii][jj] = _spikes[correct_trial_idx]
        end
    end
    spikes, poster_labels
end

function get_spatial_data(::Type{T}, celldirs::Vector{String};kwargs...) where T <: Union{SpikeCountPerGazeBin, SpikeCountPerSpatialBin}
    # TODO: This assumes that every session uses the same poster positions
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    spike_counts = Vector{Float64}[]
    trajectories = Vector{Int64}[]
    trialidx = Vector{Int64}[]
    prog = Progress(length(celldirs))
    cellcount = 0
    for (ii,sessiondir) in enumerate(sessiondirs)
        # maybe we no longer need do to this
        #jocc,qdata = cd(sessiondir) do
        #    jocc = JointOccupancy(;redo=fname->false, do_save=false, nrefinements=(p=1, g=1),trial_start=2,min_speed=1.0)
        #    qdata = UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=fname->false)
        #    jocc,qdata
        #end
        #jocc_filtered = JointFilteredOccupancy(jocc, qdata;min_place_duration=0.05, min_view_duration=0.02, min_speed=-1.0,min_place_obs=5, min_view_obs=5);
        cidx = allsessiondirs.==sessiondir
        for (jj,celldir) in enumerate(celldirs[cidx])
            obj = cd(celldir) do 
                get_num_spikes_per_bin(T;nrefinements=(p=0,g=0), trial_start=2, min_speed=1.0, min_place_duration=0.05, min_view_duration=0.02, min_place_obs=5, min_view_obs=5,kwargs...)
            end
            spike_counts_flat =  reduce(vcat, filter(sc->length(sc)>0, obj.spikecounts));
            trajectories_flat = reduce(vcat, filter(sc->length(sc)>0, obj.bins));
            _trialidx = [fill(i,length(traj)) for (i,traj) in enumerate(obj.bins)]
            trialidx_flat = reduce(vcat, filter(sc->length(sc)>0, _trialidx))
            push!(spike_counts, spike_counts_flat)
            push!(trajectories, trajectories_flat)
            push!(trialidx, trialidx_flat)
            cellcount += 1
            next!(prog)
        end
    end
    spike_counts, trajectories, trialidx
end



"""
    format_data(spikes, poster_labels)

Format the data so that spikes are encoded as number of spikes per cell per trial 
"""
function format_data(spikes, poster_labels)
    nsessions = length(spikes)
    length(poster_labels) == nsessions || error("Spikes and poster_labels should contain the same number of sessions")
    ncells = sum(length.(spikes))
    nt = maximum(length.(poster_labels))
    spikecounts = zeros(nt, ncells)
    triallabels = Vector{Vector{eltype(poster_labels[1])}}(undef, ncells)
    celloffset = 0
    for (session_spikes, session_labels) in zip(spikes, poster_labels)
        for i in 1:length(session_spikes)
            spikecounts[1:length(session_spikes[i]),celloffset+i] .= length.(session_spikes[i])
            triallabels[celloffset+i] = session_labels
        end
        celloffset += length(session_spikes)
    end
    spikecounts, triallabels
end

function format_data(spikecounts::Vector{Vector{T}}) where T <: Real
    ntrials = maximum(length.(spikecounts))
    ncells = length(spikecounts)
    scounts = zeros(ntrials,ncells)
    for (ii,_sc) in enumerate(spikecounts)
        scounts[1:length(_sc),ii] = _sc
    end
    scounts
end

function filter_min_trials(triallabel::Vector{Vector{T}};min_ntrials=10) where T
    nn = length(triallabel)
    fidx = Vector{Vector{Int64}}(undef, nn)
    for (i,label) in enumerate(triallabel)
        cc = countmap(label)
        qq = findall(cq->cq>min_ntrials,cc)
        fidx[i] = findall(in(qq), label)
    end
    fidx
end

function filter_spikecounts(spikecounts::Matrix{T}, fidx::Vector{Vector{Int64}}) where T <: Real
    ncells = size(spikecounts,2)
    nt = maximum(length.(fidx))
    spikecounts2 = zeros(T, nt, ncells) 
    for i in 1:ncells
        spikecounts2[1:length(fidx[i]),i] .= spikecounts[fidx[i],i]
    end
    spikecounts2
end

function stabilize(X::Matrix{T};dims=1) where T <: Real
    Xs = sqrt.(X)
    Xs .-= mean(Xs,dims=dims)
    Xs
end

function decode(spikecounts::Matrix{T}, triallabel::Vector{T2};ntrain=1500,ntest=100) where T2 <: Vector{T3} where T3 where T <: Real
    X = repeat(permutedims(spikecounts),1,1,1)
    Y, trainlabel, testlabel = CrossTemporalDecoding.sample_trials(X, triallabel;ntrain=ntrain,ntest=ntest)

    # convert to integer labels since this is what MultivariateStats expects
    labels = [trainlabel;testlabel]
    ulabels = unique(labels)
    idxs = Dict(l=>i for (i,l) in enumerate(ulabels))
    traintestlabel = [idxs[l] for l in labels]

    # make sure we can convert back
    @assert ulabels[traintestlabel[1:1500]] == trainlabel 
    @assert ulabels[traintestlabel[1501:end]] == testlabel

    lda = fit(SubspaceLDA, permutedims(Y[1,1:1500,:]), traintestlabel[1:1500])
    zp = predict(lda, permutedims(Y[1,1:1500,:]))
    lidx = decode(Y[1,1501:end,:], lda)
    cidx = lidx.==traintestlabel[1501:end]
    perf = sum(cidx)/length(cidx)
    decoded_labels = ulabels[lidx]
    perf, testlabel, decoded_labels, (zp, trainlabel), lda, ulabels
end

function decode(Y::Matrix{T}, lda) where T <: Real
    z = predict(lda, permutedims(Y))
    czmeans = predict(lda, lda.cmeans)
    lidx = [argmin(dropdims(sum(abs2, z[:,i] .- czmeans,dims=1),dims=1)) for i in 1:size(z,2)]
    # TODO: Convert back to origin label space
    lidx
end

function train_navigation_test_cue(spike_counts_nav, triallabels_nav, spike_counts_cue, triallabels_cue,mm::SimpleMesh;decode_goal=true, ntrain=1500,ntest=500)
    # train a spatial decoder do decode space during navgation
    X_nav = repeat(permutedims(stabilize(spike_counts_nav)),1,1,1)
    # create pseudo-population response
    Y, trainlabel, testlabel = CrossTemporalDecoding.sample_trials(X_nav, triallabels_nav;ntrain=ntrain,ntest=ntest)
    labels = [trainlabel;testlabel]
    ulabels = unique(labels)
    idxs = Dict(l=>i for (i,l) in enumerate(ulabels))
    traintestlabel = [idxs[l] for l in labels]

    lda = fit(SubspaceLDA, permutedims(Y[1,1:ntrain,:]), traintestlabel[1:ntrain])

    # create pseudo-population for cue period
    X_cue = repeat(permutedims(stabilize(spike_counts_cue)),1,1,1)
    Y_cue, trainlabel_cue, testlabel_cue = CrossTemporalDecoding.sample_trials(X_nav, triallabels_cue;ntrain=ntrain,ntest=ntest)
    # find the space spanning the variance
    pca = fit(PCA, permutedims(Y_cue[1,1:ntrain,:]))

    # find a matrix to transform the pca space in the cue period to the pca space during the navigation period
    R = procrustes_transform(permutedims(pca.proj), permutedims(lda.projw))
    # project the cue responses onto the lda space during the navigation period
    W_cue = lda.projLDA'*R*pca.proj'
    Zp = W_cue*Y_cue[1,1:ntrain,:]'
    W_nav = lda.projLDA'*lda.projw'
    #decode position
    czmeans = predict(lda, lda.cmeans)
    lidx = [argmin(dropdims(sum(abs2, Zp[:,i] .- czmeans,dims=1),dims=1)) for i in 1:size(Zp,2)] 
    decoded_labels = ulabels[lidx]  

    # find the poster positions on the mesh 
    if decode_goal
        kt = 2
    else
        kt = 1
    end
    poster_pos_on_mesh = [get_poster_position(poster_pos[poster_names[k[kt]]],mm) for k in trainlabel_cue[1:ntrain]]

    unique_poster_positions = [get_poster_position(poster_pos[poster_names[k]],mm) for k in 1:6]
    # find the performance for each poster
    perf = zeros(length(poster_pos))
    nn = fill(0, length(poster_pos))
    confusion_matrix = zeros(nelements(mm), nelements(mm))
    for (dl, tl) in zip(decoded_labels, poster_pos_on_mesh)
        ii = findfirst(unique_poster_positions.==tl)
        confusion_matrix[dl,tl] += 1
        perf[ii] += dl==tl
        nn[ii] += 1
    end
    perf./nn, confusion_matrix, unique_poster_positions, poster_names, W_nav, W_cue
end

function get_colors(traj::Vector{Tuple{Int64, Int64}})
    nn = length(traj)
    colors = Vector{RGB{Colors.FixedPointNumbers.N0f8}}(undef, nn)
    for (ii,(p1,p2)) in enumerate(traj)
        c1 = parse(Colorant, poster_color[poster_names[p1]])
        c2 = parse(Colorant, poster_color[poster_names[p2]])
        c12 = weighted_color_mean(0.7, c1, c2)
        colors[ii] = c12
    end
    colors
end

function get_confusion_matrix(true_labels::AbstractVector{T}, decoded_labels::AbstractVector{T}) where T <: Integer
    lmax = max(maximum(true_labels), maximum(decoded_labels))
    C = zeros(lmax, lmax)
    for (tl,dl) in zip(true_labels, decoded_labels)
        C[dl,tl] += 1.0
    end
    C
end


## plots
function plot_decoder_space(zq::Matrix{<:Real}, triallabels::Vector{Tuple{Int64, Int64}};_plot_theme=plot_theme)
    tcolor = get_colors(triallabels) 
    unique_labels = unique(triallabels)
    # hard coded; maybe make this a bit more gneeral
    img = Matrix{eltype(tcolor)}(undef, 6, 6)
    for ul in unique_labels
        ii = findfirst(tt->tt==ul, triallabels)
        img[ul[1], ul[2]] = tcolor[ii]
    end
    w = 1024
    h = (9/16)*w
    with_theme(_plot_theme) do
        fig = Figure(size=(w,h))
        ax = Axis3(fig[1,1])
        ax.xticklabelsvisible = false
        ax.yticklabelsvisible = false
        ax.zticklabelsvisible = false
        ax.xlabelvisible = false
        ax.ylabelvisible = false
        ax.zlabelvisible = false
        scatter!(ax, Point3f.(eachcol(zq)),color=tcolor)
        bbox = lift(fig.scene.viewport) do vp
            w,h = vp.widths
            x0,y0 = vp.origin
            BBox(w-250, w-50, h-250, h-50)
        end
        ax2 = Axis(fig, bbox=bbox)
        image!(ax2, img;interpolate=false)
        ax2.xticks = ([0:5;] .+ 0.5, String.(poster_names))
        ax2.yticks = ([0:5;] .+ 0.5, String.(poster_names))
        ax2.xticklabelsize = 16
        ax2.yticklabelsize = 16
        ax2.xticklabelrotation = -π/6
        fig
    end
end

function plot_spatial_performance(true_labels::Vector{T}, decoded_labels::Vector{T},mm::SimpleMesh;_plot_theme=plot_theme, kwargs...) where T <: Integer
    Z = zeros(nelements(mm)) 
    nn = fill(0, nelements(mm))
    colormap = get(kwargs, :colormap, :navia)
    for (tl,dl) in zip(true_labels, decoded_labels)
        Z[tl] += tl==dl
        nn[tl] += 1
    end
    Z ./= nn
    with_theme(_plot_theme) do
        fig = Figure()
        ax = Axis(fig[1,1], aspect=1.0)
        hidedecorations!(ax)
        ax.bottomspinevisible = false
        ax.leftspinevisible = false
        viz!(ax, mm;color=Z,colormap=colormap)
        text!(ax, Point2f.(poster_pos[k] for k in keys(poster_pos)), text=String.(collect(keys(poster_pos))), 
                        align=(:center, :baseline), color=:gray85)
        Colorbar(fig[1,2], colorrange=extrema(filter(isfinite, Z)), label="Performance", colormap=colormap)
        fig
    end
end

function plot_spatial_confusion_matrix(confusion_matrix, mm::SimpleMesh;_plot_theme=plot_theme,colormap=:managua)
    perf = diag(confusion_matrix)./dropdims(sum(confusion_matrix,dims=1),dims=1)
    with_theme(_plot_theme) do 
        fig = Figure()
        lscene = LScene(fig[1,1], show_axis=false)
        mm_upper = Translate(0.0,0.0, 10.0)(mm)
        mm_lower = Translate(0.0,0.0, -10.0)(mm)
        viz!(lscene, mm_upper;showsegments=true, color=:gray20)
        viz!(lscene, mm_lower;showsegments=true,color=perf, colormap=colormap)
        ls = Tuple{Point3f, Point3f}[]
        lw = Float64[]
        cc = Float64[]
        for j in 1:size(confusion_matrix,2)
            p2 = Point3f(Tuple(centroid(mm_lower[j])))
            for i in 1:size(confusion_matrix,1)
                if confusion_matrix[i,j] > 0
                    p1 = Point3f(Tuple(centroid(mm_upper[i])))
                    push!(ls,(p2,p1) )
                    push!(lw, confusion_matrix[i,j])
                    push!(cc, perf[j])
                end
            end
        end
        lw = 5*(lw .- minimum(lw))./(maximum(lw) - minimum(lw)) .+ 0.1
        linesegments!(lscene, ls, linewidth=lw, color=cc, colormap=colormap)
        text!(lscene, Point3f.((poster_pos[k]...,-9) for k in keys(poster_pos)), text=String.(collect(keys(poster_pos))), 
                        align=(:center, :baseline), color=:salmon, overdraw=true)#, markerspace=:data, fontsize=1.5)
        linesegments!(lscene, [(Point3f(poster_pos[k]..., -10), Point3f(poster_pos[k]...,-9)) for k in keys(poster_pos)], color=:salmon,
                            overdraw=true)
        Colorbar(fig[1,2], colorrange=extrema(filter(isfinite, perf)), colormap=colormap, label="Performance")
        fig
    end
end

function plot_decoder_contribution(W_nav::Matrix{T}, W_cue::Matrix{T};_plot_theme = plot_theme, cell_color::Union{AbstractVector{T2}, Nothing}=nothing, color_legend::Dict{Symbol, Symbol}=Dict{Symbol, Symbol}()) where T <: Real where T2
    if cell_color === nothing
        cell_color = :royalblue2
    end
    w_nav = dropdims(sum(abs2, W_nav,dims=1),dims=1)
    w_cue = dropdims(sum(abs2, W_cue,dims=1),dims=1)
    with_theme(_plot_theme) do
        fig = Figure()
        ax = Axis(fig[1,1])
        scatter!(ax, w_nav, w_cue, color=cell_color)
        ax.xlabel = "Code weight during navigation"
        ax.ylabel = "Code weight during cue"
        if !isempty(color_legend)
            Legend(fig[1,1], [MarkerElement(marker=:cicle, color=color_legend[k]) for k in keys(color_legend)],
                            String.(collect(keys(color_legend))),valign=:top, halign=:right, tellwidth=false,
                            tellheight=false)
        end
        fig
    end
end