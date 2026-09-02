"""
Look at a cells response when the gaze is in its view field overlapping with a goal poster
"""
function analyze_goal_poster_selectivity(;kwargs...)
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
        return (X_goal, w_goal), (X_nongoal, w_nongoal), poster_selectivity
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
    (X_goal, w_goal), (X_nongoal, w_nongoal), poster_selectivity
end

## plots

function plot_goal_vs_non_goal_poster_tuning(X_goal, w_goal, X_nongoal, w_nongoal)
    λ_goal = X_goal./w_goal
    fidx_goal = [findall(w_goal[:,i] .> 0.02) for i in 1:size(w_goal,2)]
    λ_nongoal = X_nongoal./w_nongoal
    fidx_nongoal = [findall(w_nongoal[:,i] .> 0.02) for i in 1:size(w_nongoal,2)]

    @show fidx_goal    
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
    @show size(μ_nongoal)
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