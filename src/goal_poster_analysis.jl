struct GoalPosterSelectivity{T<:GazeResponseFieldsAll}
    rf_gaze::T
    X_goal::Matrix{Float64}
    w_goal::Matrix{Float64}
    X_nongoal::Matrix{Float64}
    w_nongoal::Matrix{Float64}
    poster_selectivity::Matrix{Bool}
    args::Dict{Symbol,Any}
end

DPHT.filename(::Type{GoalPosterSelectivity}) = "goal_poster_selectivity.jld2"

function process_kwargs(::Type{GoalPosterSelectivity{T}},h::UInt32=zero(UInt32);exclude_origin_trials=false, kwargs...)  where T <: GazeResponseFieldsAll
    h = process_kwargs(T, h;kwargs...)
    if exclude_origin_trials
        h = CRC32c.crc32c(string(:exclude_origin_trials=>exclude_origin_trials),h)
    end
    h
end

function get_rate(X::GoalPosterSelectivity)
    X.X_goal./X.w_goal, X.X_nongoal./X.w_nongoal
end

function generate_surrogates(X::GoalPosterSelectivity;nruns=1000)
    λ_goal = X.X_goal./X.w_goal
    statstest = fill(false, size(λ_goal,2))
    fidx_goal = [findall(X.w_goal[:,i] .> 0.02) for i in 1:size(X.w_goal,2)]
    if isempty(fidx_goal)
        return statstest
    end
    λ_nongoal = X.X_nongoal./X.w_nongoal
    fidx_nongoal = [findall(X.w_nongoal[:,i] .> 0.02) for i in 1:size(X.w_nongoal,2)]

    μ_goal = [mean(λ_goal[fidx_goal[i],i]) for i in 1:size(X.w_goal,2)]
    μ_nongoal = [mean(λ_nongoal[fidx_nongoal[i],i]) for i in 1:size(X.w_nongoal,2)]
    qidx = findall(isfinite,μ_goal)
    statstest = fill(false, length(μ_goal))

    Δ = μ_goal .- μ_nongoal 
    Δs = zeros(1000, length(Δ))
    for i in axes(μ_nongoal,1)
        _fidx_g = fidx_goal[i]
        _fidx_ng = fidx_nongoal[i]
        nn_g = length(_fidx_g)
        nn_ng = length(_fidx_ng)
        if nn_g == 0 || nn_ng == 0
            continue
        end
        # grab random indices
        for j in 1:1000 
            # grab random nn_g 
            _μg = mean([rand() < 0.5 ? λ_goal[rand(_fidx_g),i] : λ_nongoal[rand(_fidx_ng),i] for _ in 1:nn_g])
            _μng = mean([rand() < 0.5 ? λ_goal[rand(_fidx_g),i] : λ_nongoal[rand(_fidx_ng),i] for _ in 1:nn_ng])
            Δs[j,i] = _μg - _μng
        end
    end
    Δ, Δs
end
function Hippocampus.issignificant(X::GoalPosterSelectivity;pv_threshold=0.05, side=:both)
    λ_goal = X.X_goal./X.w_goal
    statstest = fill(false, size(λ_goal,2))
    fidx_goal = [findall(X.w_goal[:,i] .> 0.02) for i in 1:size(X.w_goal,2)]
    if isempty(fidx_goal)
        return statstest
    end
    λ_nongoal = X.X_nongoal./X.w_nongoal
    fidx_nongoal = [findall(X.w_nongoal[:,i] .> 0.02) for i in 1:size(X.w_nongoal,2)]

    μ_goal = [mean(λ_goal[fidx_goal[i],i]) for i in 1:size(X.w_goal,2)]
    μ_nongoal = [mean(λ_nongoal[fidx_nongoal[i],i]) for i in 1:size(X.w_nongoal,2)]
    qidx = findall(isfinite,μ_goal)
    statstest = fill(false, length(μ_goal))

    Δ = μ_goal .- μ_nongoal 
    Δs = zeros(1000, length(Δ))
    for i in axes(μ_nongoal,1)
        _fidx_g = fidx_goal[i]
        _fidx_ng = fidx_nongoal[i]
        nn_g = length(_fidx_g)
        nn_ng = length(_fidx_ng)
        if nn_g == 0 || nn_ng == 0
            continue
        end
        # grab random indices
        for j in 1:1000 
            # grab random nn_g 
            _μg = mean([rand() < 0.5 ? λ_goal[rand(_fidx_g),i] : λ_nongoal[rand(_fidx_ng),i] for _ in 1:nn_g])
            _μng = mean([rand() < 0.5 ? λ_goal[rand(_fidx_g),i] : λ_nongoal[rand(_fidx_ng),i] for _ in 1:nn_ng])
            Δs[j,i] = _μg - _μng
        end
        Δu = percentile(Δs[:,i],100*(1-pv_threshold))
        Δl = percentile(Δs[:,i],100*(pv_threshold))
        if side == :both
            statstest[i] = (Δ[i] > Δu) || (Δ[i] < Δl)
        elseif side == :left
            statstest[i] = Δ[i] < Δl
        else
            statstest[i] = Δ[i] > Δu
        end
    end
    return statstest
end
"""
Look at a cells response when the gaze is in its view field overlapping with a goal poster
"""
function analyze_goal_poster_selectivity(;exclude_origin_trials=false, kwargs...)
    nrefinments = get(kwargs, :nrefinements, (p=3,g=2))
    mm = get_maze_mesh(;nrefinements=nrefinments.g)
    rf_gaze = get_response_fields(GazeResponseFields, get(kwargs, :nshuffles, 1000);kwargs...)
    jm = JointMap(;kwargs...)
    udata = cd(DPHT.process_level("session")) do
        UnityData()
    end
    viewfields = getfields(rf_gaze;cluster_threshold=0.001)
    # find the overlap between the view fields and the poster
    poster_selectivity = fill(false, 6, length(viewfields))
    for (ii,vidx) in  enumerate(viewfields)
        poster_selectivity[:,ii] = find_poster_view_intersection(vidx, mm)
    end
    nt = numtrials(udata)
    X_goal = zeros(nt, length(viewfields))
    w_goal = zeros(nt, length(viewfields))
    X_nongoal = zeros(nt, length(viewfields))
    w_nongoal = zeros(nt, length(viewfields))
    if !(any(poster_selectivity))
        return (X_goal, w_goal), (X_nongoal, w_nongoal), poster_selectivity, rf_gaze
    end
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        if w == 0
            continue
        end
        vidx,pidx,hidx,tidx = Tuple(qidx)
        # are we in the view field?
        kk = findfirst([in(vf)(vidx) for vf in viewfields])
        if kk === nothing
            continue
        end
        
        # is there a poster in this field
        poster_pref = findall(poster_selectivity[:,kk])
        if isempty(poster_pref)
            continue
        end
        # check if this is a correct trial
        if 30 < udata.triggers[tidx,3] < 40
            if tidx > 1 && exclude_origin_trials
                # we do not want to include trials for which the previous trial used the
                # poster_pref as a target, since the subject will always be moving away from the
                # poster in the current trial and we do not want that to influence the no-goal firing rate
                if in(poster_pref)(udata.triggers[tidx-1,1] - 10)
                    continue
                end
            end
            posterid = udata.triggers[tidx,1] - 10
            if in(poster_pref)(posterid)
                X_goal[tidx,kk] += w
                w_goal[tidx,kk] += occ
            else
                X_nongoal[tidx, kk] +=  w
                w_nongoal[tidx,kk] += occ 
            end
        end
    end
    (X_goal, w_goal), (X_nongoal, w_nongoal), poster_selectivity, rf_gaze
end

function GoalPosterSelectivity(::Type{T};do_save=true, redo=fanme->false, kwargs...) where T <: GazeResponseFieldsAll
    fname = DPHT.filename(GoalPosterSelectivity)
    h = process_kwargs(GoalPosterSelectivity{T};kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = false
    if !isfile(fname) || redo(fname)
        do_compute=true
    end
    if do_compute
        (X_goal, w_goal), (X_non_goal, w_nongoal), poster_selectivity,rf_gaze = analyze_goal_poster_selectivity(;kwargs...)
        dd = Dict(kwargs)
        obj = GoalPosterSelectivity{T}(rf_gaze,X_goal, w_goal, X_non_goal, w_nongoal, poster_selectivity, dd)
        if do_save
            save_jld2(obj, fname)
        end
    else
        obj = load_jld2(GoalPosterSelectivity{T}, fname)
    end
    obj
end

## plots

function plot_goal_poster_selectivity!(lg, gs::GoalPosterSelectivity;label=["A","B"],kwargs...)
    λ_goal = gs.X_goal./gs.w_goal
    fidx_goal = [findall(gs.w_goal[:,i] .> 0.02) for i in 1:size(gs.w_goal,2)]
    μ_goal = [mean(λ_goal[fidx_goal[i],i]) for i in 1:size(gs.w_goal,2)]
    λ_nongoal = gs.X_nongoal./gs.w_nongoal
    fidx_nongoal = [findall(gs.w_nongoal[:,i] .> 0.02) for i in 1:size(gs.w_nongoal,2)]
    μ_nongoal = [mean(λ_nongoal[fidx_nongoal[i],i]) for i in 1:size(gs.w_goal,2)]

    Δ, Δs = generate_surrogates(gs)
    ax0 = Axis(lg[1,1])
    Label(lg[1,1,TopLeft()], label[1])
    # For each view view, plot the goal and no-goal firing rate next to each other
    ddp = vec(permutedims([fill(1, length(μ_goal)) fill(2, length(μ_nongoal))]))
    yyp = vec(permutedims([μ_goal μ_nongoal]))
    xxp = vec(permutedims([1:length(μ_goal) 1:length(μ_nongoal)]))
    xt = vec(permutedims([1 2] .+ 2*([1:length(μ_goal);] .-1)))
    barplot!(ax0, xt, yyp, direction=:x)
    ax0.yticks = (xt, repeat(["Goal", "No goal"], length(μ_goal)))
    ax0.yticklabelrotation = -π/6
    ax0.yticklabelalign = (:right, :center)
    ax0.xlabel = "Firing rate [Hz]"
    ax0.yticksvisible = true

    ax = Axis(lg[1,2])
    Label(lg[1,2, TopLeft()], label[2])
    yy = vec(Δs[:,isfinite.(μ_goal)])
    xx = reduce(vcat, [fill(i,size(Δs,1)) for i in findall(isfinite.(μ_goal))])
    boxplot!(ax, xx, yy;show_outliers=false, show_notch=true, orientation=:horizontal)
    scatter!(ax, Δ[isfinite.(μ_goal)], [1:sum(isfinite.(μ_goal));] ;color=:orange)
    ax.yticklabelsvisible = false
    ax.yticksvisible = false
    ax.xlabel = "Rate(goal) - Rate(no-goal)"
    for _ax in [ax,ax0]
        _ax.xticksvisible = true
        _ax.xticklabelsvisible = get(kwargs, :xticklabelsvisible, true)
        _ax.xlabelvisible = get(kwargs, :xlabelvisible, true)
    end
    # add an axis showing the poster labels
    ax1 = Axis(lg[1,3])
    imgs = [load(poster_img[k]) for k in poster_names]
    posterpref = findfirst.(eachcol(gs.poster_selectivity))
    qidx = posterpref.!== nothing
    scatter!(ax1, fill(1.0,sum(qidx)), findall(qidx), marker=imgs[posterpref[qidx]], markersize=45)
    linkyaxes!(ax, ax1)
    hidedecorations!(ax1)
    hidespines!(ax1)
    colsize!(lg, 3, 50)
end

function plot_goal_poster_selectivity(gs::GoalPosterSelectivity;kwargs...)
    fidx_goal = [findall(gs.w_goal[:,i] .> 0.02) for i in 1:size(gs.w_goal,2)]
    nn = sum(length.(fidx_goal).>0)
    width = 600
    with_theme(theme_minimal()) do
        fig = Figure(size=(width,250))
        lg = GridLayout(fig[1,1])
        plot_goal_poster_selectivity!(lg, gs;kwargs...)
        fig
    end
end

function plot_goal_vs_non_goal_poster_tuning(X_goal, w_goal, X_nongoal, w_nongoal)
    λ_goal = X_goal./w_goal
    fidx_goal = [findall(w_goal[:,i] .> 0.02) for i in 1:size(w_goal,2)]
    λ_nongoal = X_nongoal./w_nongoal
    fidx_nongoal = [findall(w_nongoal[:,i] .> 0.02) for i in 1:size(w_nongoal,2)]

    μ_goal = [mean(λ_goal[fidx_goal[i],i]) for i in 1:size(w_goal,2)]

    μ_nongoal = zeros(1000, length(μ_goal))
    qidx = findall(isfinite,μ_goal)

    statstest = fill(false, length(μ_goal))
    for i in axes(μ_nongoal,2)
        μ_nongoal[:,i] = [mean(λ_nongoal[shuffle(fidx_nongoal[i])[1:length(fidx_goal[i])],i]) for _ in 1:1000]
        if isfinite(μ_goal[i])
            statstest[i] = μ_goal[i] > percentile(μ_nongoal[:,i], 85)
        end
    end
    @show statstest
    with_theme(theme_minimal()) do
        fig = Figure(size=(100,300))
        ax = Axis(fig[1,1])
        yy =vec(μ_nongoal[:,qidx])
        xx = vec([fill(1.0, 1000) fill(2.0, 1000)][:,qidx])
        boxplot!(ax, xx, yy,show_outliers=false, show_notch=true)
        scatter!(ax, qidx, μ_goal[qidx], color=:orange)
        ax.xticksvisible = false
        ax.xticklabelsvisible = false
        ax.ylabel = "Firing rate [Hz]"
        fig
    end
end