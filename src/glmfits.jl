using SpecialFunctions
using Random
using Optim
using GLM
using ReverseDiff
using ADTypes

struct GLMFit
    β_pos::Vector{Float64}
    deviance_pos::Tuple{Float64, Float64}
    β_gaze::Vector{Float64}
    deviance_gaze::Tuple{Float64, Float64}
    β_hd::Vector{Float64}
    deviance_hd::Tuple{Float64, Float64}
    β_all::Vector{Float64}
    deviance_all::Tuple{Float64, Float64}
end

DPHT.filename(::Type{GLMFit}) = "glmfit.jld2"
DPHT.level(::Type{GLMFit}) = "cell"

struct GLMFitH{N}
    β::Array{Float64,3}
    ll::Matrix{Float64}
    α::Vector{Float64}
    trainidx::Matrix{Int64}
    nspikes::Vector{Int16}
    dims::NTuple{N,Symbol}
    qidx::Vector{CartesianIndex{4}}
    converged::Bool
end

function GLMFitH{N}(β, ll, α, trainidx, nspikes, dims::NTuple{N,Symbol}, qidx) where N
    GLMFitH{N}(β, ll, α, trainidx, nspikes, dims, qidx, true)
end

function GLMFitH(β, ll, α, trainidx, nspikes, dims::NTuple{N,Symbol}, qidx) where N
    GLMFitH{N}(β, ll, α, trainidx, nspikes, dims, qidx, true)
end

DPHT.filename(::Type{GLMFitH{N}}) where N = "glmfith.jld2"

function DPHT.filename(::Type{GLMFitH{N}},dims::NTuple{N,Symbol}) where N
    nn = join(String.(dims),'_')
    fname = "glmfith_$(nn).jld2"
    fname
end

DPHT.level(::Type{GLMFitH{N}}) where N = "cell"

function process_kwargs(::Type{GLMFitH{N}};α=10.0.^[-2,-3,-4,-5,-6], nruns=10, kwargs...) where N
    h = UInt32(0)
    h = CRC32c.crc32c(string((:α=>α)),h)
    h = CRC32c.crc32c(string((:nruns=>nruns)),h)
    h
end

function plot_glmfit(glmfit::GLMFitH{N},args...;kwargs...) where N
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_glmfit!(lg, glmfit,args...;kwargs...)
        fig
    end
end

function plot_glmfit!(lg, glmfit_joint::GLMFitH{N}, glmfit::Vector{GLMFitH{1}}) where N
    nspikes = glmfit_joint.nspikes
    nruns = size(glmfit_joint.ll,1)
    trainidx = glmfit_joint.trainidx
    ntest = length(nspikes) - size(trainidx,1)
    ll = zeros(nruns)
    for (i,_glmfit) in enumerate(glmfit)
        aidx = findfirst(_glmfit.α.==glmfit_joint.α[i])
        β = _glmfit.β[:,:,aidx]
        if _glmfit.dims[1] == :p
            mm = floor_topology3()
            nb = nelements(mm)
            m = 2
        elseif _glmfit.dims[1] == :g
            mm = get_maze_mesh()
            nb = nelements(mm)
            m = 1
        else
            nb = 24
            m = 3
        end
        X = zeros(nb+1, ntest)
        X[nb+1,:] .= 1.0
        for (j,_trainidx) in enumerate(eachcol(trainidx))
            testidx = setdiff(1:length(nspikes), _trainidx)
            for (k,qq) in enumerate(glmfit_joint.qidx[testidx])
                X[qq.I[m],k] = 1.0
            end
            ll[j] += logprob(β[:,j], X, nspikes[testidx])
        end
    end

    ax = Axis(lg[1,1])
    scatter!(ax, glmfit_joint.ll[:,1], ll)
    ablines!(ax, 0.0, 1.0, linestyle=:dot, color=:black)
    ax.xlabel = "log-likelihood joint"
    ax.ylabel = "log-likelihood ind"
end

function plot_glmfit!(lg, glmfit::GLMFitH{N},aidx::Union{Int64, Nothing}=nothing;kwargs...) where N
    #show log-likelihoods for each cross-validation fold and the resulting map
    if aidx === nothing
        α,aidx = get_best_α(glmfit)
    end
    if :g in(glmfit.dims)
        pax = LScene 
        mm = get_maze_mesh()
    else
        pax = Axis
        mm = floor_topology3()
    end
    ll = glmfit.ll
    nspikes = glmfit.nspikes
    # establish null likelihood
    ll0 = zeros(size(ll,1))
    for i in 1:length(ll0)
        trainidx = glmfit.trainidx[:,i]
        testidx = setdiff(1:length(nspikes), trainidx)
        P = fit(Poisson,nspikes[trainidx])
        ll0[i] = mean(logpdf.(P, nspikes[testidx]))
    end
    xx = repeat(glmfit.α, 1,size(ll,1))[:]
    yy = permutedims(glmfit.ll)[:]
    points = Point2f[]
    for i in 1:size(ll,1)
        for j in 1:size(ll,2)
            push!(points, Point2f(glmfit.α[j],ll[i,j]))
        end
        push!(points, Point2f(NaN))
    end
    lg2 = lg[1,1]
    ax1 = Axis(lg2[1,1],xscale=log10)
    Label(lg2[1,1,TopLeft()],"A")
    lines!(ax1, points)
    scatter!(ax1, xx,yy)
    if get(kwargs, :ylabelvisible, true)
        ax1.ylabel = "Mean log-likelihood"
    end
    if get(kwargs, :xlabelvisible, true)
        ax1.xlabel = "Smoothing factor"
    end
    vlines!(ax1, glmfit.α[aidx], color=:black, linestyle=:dot)
    ax2 = Axis(lg2[2,1])
    Label(lg2[2,1,TopLeft()],"B")
    scatter!(ax2, ll[:,aidx], ll0)
    ax2.xlabel = "Log-likelihood"
    ax2.ylabel = "Log-likelihood null"
    ablines!(ax2, 0.0, 1.0, linestyle=:dot, color=:black)
    ax2.xticklabelrotation = -π/7
    axp = pax(lg[1,2])
    Label(lg[1,2,TopLeft()],"C")
    if isa(axp, LScene)
        plotmesh!(axp, mm;showsegments=false, color=color=glmfit.β[1:end-1, 1,aidx], ceiling_offset=10, floor_offset=-20,shading=false)
    elseif first(glmfit.dims) == :hd
        θ = range(0.0, stop=2π, length=24)
        lines!(axp, θ, glmfit.β[1:end-1,1,aidx])
    else
        viz!(axp, mm;showsegments=false, color=glmfit.β[1:end-1,1,aidx])
    end
    ax1,axp
end

function get_best_α(glmfit::GLMFitH{N}) where N
    nruns,nα = size(glmfit.ll)
    nn = zeros(Int64, nα)
    for i in 1:nα
        for j in 1:nα
            if i == j
                continue
            end
            nn[i] += sum(glmfit.ll[:,i] .> glmfit.ll[:,j])
        end
    end
    idx = argmax(nn)
    glmfit.α[idx],idx
end

"""
Get the probability that a null model produces the log-likelihoods found in 
cross-validation fold `idx`.
"""
function get_significance_level(glmfit::GLMFitH{N},idx::Integer;α=0.05) where N
    ll = glmfit.ll[:,idx]
    nspikes = glmfit.nspikes
    ll0 = zeros(length(ll))
    for i in 1:length(ll0)
        trainidx = glmfit.trainidx[:,i]
        testidx = setdiff(1:length(nspikes), trainidx)
        P = fit(Poisson,nspikes[trainidx])
        ll0[i] = mean(logpdf.(P, nspikes[testidx]))
    end
    nn = sum(ll .> ll0)
    pv = 1 - cdf(Binomial(length(ll),α),nn)
    nn, pv
end

function get_significance_level(glmfit::GLMFitH{N};α=0.05) where N
    nruns = size(glmfit.ll,1)
    trainidx = glmfit.trainidx
    if length(glmfit.dims) > 1
        # compare signfiance level to the product of the individuals
        ll = zeros(nruns) 
        for (i,d) in enumerate(glmfit.dims)
            # need to load the corresponding object
            files = glob("glmfith_$(d)_*.jld2")
            for f in files
                _glmfit = JLD2.load(f)
                if size(_glmfit["data"].ll,1) == nruns
                    aidx = findfirst(_glmfit["data"].α .== glmfit.α[i])
                    if aidx !== nothing
                        β = _glmfit["data"].β[:,:,aidx]

                        ll .+= _glmfit["data"].ll[:,aidx]

                        break
                    end
                end
            end
        end
        nn = sum(glmfit.ll[:,1] .> ll)
        pv = 1 - cdf(Binomial(nruns,α), nn)
        return nn, pv
    else
        aidx = get_best_α(glmfit)
        get_significance_level(glmfit, aidx;α=α)
    end

end

function get_num_spikes(vpvrp::ViewAndPlaceRepresentationNew, vpoc::ViewAndPlaceOccupancy)
    mm = get_maze_mesh()
    m_floor = floor_topology3()
    # use 24 bins for head direction
    hd_bins = range(0.0, stop=2π, length=24)
    nt = size(vpoc.weight_view,3)
    cc = zeros(Int64, nelements(mm),nelements(m_floor),nt)
    for i in 1:size(cc,3)
        for _idx in vpvrp.placeviewidx[i]
            vidx = vpoc.viewbin_idx[i][_idx]
            pidx = vpoc.placebin_idx[i][_idx]
            #θidx = argmax(cos.(hd_bins .- ))
            if (vidx > 0) && (pidx>0)
                cc[vidx,pidx,i] += 1
            end
        end
    end
    cc
end

function get_num_spikes(vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy)
    mm = get_maze_mesh()
    m_floor = floor_topology3()
    # use 24 bins for head direction
    nt = length(vpvrp.events)
    cc = Dict{CartesianIndex{4}, Int16}()
    for i in 1:nt
        for _idx in  vpvrp.placeviewidx[i]
            aidx = jocc.index[i][_idx]
            vidx,pidx,hidx = Tuple(aidx) 
            if (vidx == 0 || pidx ==0) || (hidx==0)
                continue
            end
            qq = CartesianIndex(vidx,pidx,hidx,i)
            cc[qq] = get(cc, qq, zero(Int16)) + one(Int16)
        end
    end
    cc
end

"""
Session level
"""
function get_data(vpoc::ViewAndPlaceOccupancy, unity_gaze_data::UnityRaytraceData)
    # find place bins with minimum amount of occpuancy
    m_floor = floor_topology3()
    goodbinidx = findall(dropdims(sum(dropdims(sum(vpoc.weight_view,dims=1),dims=1) .> 0.02,dims=2),dims=2).> 5)
    # get the speed per place bin
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, vpoc.placebin_idx,nelements(m_floor))
    qidx1 = vpoc.weight_view[:,goodbinidx,:].>0
    qidx2 = reshape(vv[goodbinidx,:],1,length(goodbinidx), size(vv,2)).>1.0
    qidx = findall(qidx1.&qidx2)
    qidx2 = [CartesianIndex(ii.I[1], goodbinidx[ii.I[2]], ii.I[3]) for ii in qidx]
    np = sum(length.(unity_gaze_data.head_direction))
    pos = zeros(2,np)
    gaze = zeros(3,np)
    hd = zeros(np)
    offset = 0
    nt = size(vpoc.weight_view,3)
    flat_idx = fill(0, 2, np)
    for tidx in 1:nt 
        _qidx = qidx2[getindex.(qidx2,3) .== tidx]

        pidx = getindex.(_qidx, 2)
        vidx = getindex.(_qidx,1)

        trialpos = unity_gaze_data.position[tidx]
        trialhd = unity_gaze_data.head_direction[tidx]
        trialgaze = unity_gaze_data.gaze[tidx]

        poidx = in(pidx).(vpoc.placebin_idx[tidx])
        goidx = in(vidx).(vpoc.viewbin_idx[tidx])

        pgidx = poidx.&goidx
        npg = sum(pgidx)
        pos[:,offset+1:offset+npg] .= trialpos[1:2,pgidx]
        gaze[:,offset+1:offset+npg] .= trialgaze[:,pgidx]
        hd[offset+1:offset+npg] .= trialhd[pgidx]
        flat_idx[1,offset+1:offset+npg] .= vpoc.placebin_idx[tidx][pgidx]
        flat_idx[2,offset+1:offset+npg] .= vpoc.viewbin_idx[tidx][pgidx]
        offset += npg
    end
    pos[:,1:offset], hd[1:offset], gaze[:,1:offset], flat_idx[:,1:offset], qidx2
end

function get_data(jocc::JointOccupancy, unity_gaze_data::UnityRaytraceData)
    # find place bins with minimum amount of occpuancy
    nt = numtrials(unity_gaze_data)
    m_floor = floor_topology3()
    place_weight = zeros(nelements(m_floor),nt)
    placebin_idx = Vector{Vector{Int64}}(undef, nt)
    for (k,v) in jocc.weight
        vidx,pidx,hidx,tidx = Tuple(k)
        place_weight[pidx,tidx] += v
    end
    for i in 1:nt
        aidx = jocc.index[i]
        placebin_idx[i] = fill(0, length(aidx))
        for (j,kk) in enumerate(aidx)
            placebin_idx[i][j] = kk.I[2]
        end
    end
    goodbinidx = findall(dropdims(sum(place_weight .> 0.02,dims=2),dims=2).> 5)
    ff = in(goodbinidx)
    # get the speed per place bin
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, placebin_idx,nelements(m_floor))
    qidx = CartesianIndex{4}[]
    for k in keys(jocc.weight)
        pidx = k.I[2]
        tidx = k.I[4]
        if ff(pidx) && (vv[pidx,tidx]  > 1.0)
            push!(qidx, k)
        end
    end
    np = sum(length.(unity_gaze_data.head_direction))
    pos = zeros(2,np)
    gaze = zeros(3,np)
    hd = zeros(np)
    offset = 0
    flat_idx = Vector{CartesianIndex{3}}(undef, np)
    for tidx in 1:nt 
        _qidx = qidx[getindex.(qidx,3) .== tidx]
        _qidx = [CartesianIndex(kk.I[1], kk.I[2], kk.I[3]) for kk in _qidx]

        trialpos = unity_gaze_data.position[tidx]
        trialhd = unity_gaze_data.head_direction[tidx]
        trialgaze = unity_gaze_data.gaze[tidx]

        pgidx = in(_qidx).(jocc.index[tidx])

        npg = sum(pgidx)
        pos[:,offset+1:offset+npg] .= trialpos[1:2,pgidx]
        gaze[:,offset+1:offset+npg] .= trialgaze[:,pgidx]
        hd[offset+1:offset+npg] .= trialhd[pgidx]
        flat_idx[offset+1:offset+npg] .= jocc.index[tidx][pgidx]
        offset += npg
    end
    pos[:,1:offset], hd[1:offset], gaze[:,1:offset], flat_idx[1:offset], qidx
end

function fit_glm(vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy,unity_gaze_data::UnityRaytraceData)
    nt = numtrials(unity_gaze_data)
    t1 = time()
    m_floor = Shadow("xy")(floor_topology3())
    mm = get_maze_mesh()
    m_floor = floor_topology3()
    hdbins = range(0.0, stop=2π,length=24)
    place_weight = zeros(nelements(m_floor),nt)
    placebin_idx = Vector{Vector{Int64}}(undef, nt)
    for (k,v) in jocc.weight
        vidx,pidx,hidx,tidx = Tuple(k)
        place_weight[pidx,tidx] += v
    end
    for i in 1:nt
        aidx = jocc.index[i]
        placebin_idx[i] = fill(0, length(aidx))
        for (j,kk) in enumerate(aidx)
            placebin_idx[i][j] = kk.I[2]
        end
    end
    goodbinidx = findall(dropdims(sum(place_weight .> 0.02,dims=2),dims=2).> 5)
    ff = in(goodbinidx)
    # get the speed per place bin
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, placebin_idx,nelements(m_floor))
    qidx = CartesianIndex{4}[]
    for k in keys(jocc.weight)
        pidx = k.I[2]
        tidx = k.I[4]
        if ff(pidx) && (vv[pidx,tidx]  > 1.0)
            push!(qidx, k)
        end
    end
    # replace with goodbinidx 
    mpos = stack(Tuple.(centroid.(m_floor[getindex.(qidx,2)])))
    mgaze = stack(Tuple.(centroid.(mm[getindex.(qidx,1)])))
    mhd = hdbins[getindex.(qidx,3)] 
    t2 = time()
    cc = get_num_spikes(vpvrp, jocc)
    nspikes = zeros(Int16, length(qidx))
    for (ii,k) in enumerate(qidx)
        if k in keys(cc)
            nspikes[ii] = cc[k]
        end
    end
    nspikes, mpos, mgaze,mhd, qidx
end

function fit_glm(vpvrp::ViewAndPlaceRepresentationNew, vpoc::ViewAndPlaceOccupancy,unity_gaze_data::UnityRaytraceData)
    t1 = time()
    m_floor = Shadow("xy")(floor_topology3())
    mm = get_maze_mesh()
    goodbinidx = findall(dropdims(sum(dropdims(sum(vpoc.weight_view,dims=1),dims=1) .> 0.02,dims=2),dims=2).> 5)
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, vpoc.placebin_idx,nelements(m_floor))
    qidx1 = vpoc.weight_view[:,goodbinidx,:].>0
    qidx2 = reshape(vv[goodbinidx,:],1,length(goodbinidx), size(vv,2)).>1.0
    qidx = findall(qidx1.&qidx2)
    # replace with goodbinidx 
    qidx2 = [CartesianIndex(ii.I[1], goodbinidx[ii.I[2]], ii.I[3]) for ii in qidx]
    mpos = stack(Tuple.(centroid.(m_floor[getindex.(qidx,2)])))
    mgaze = stack(Tuple.(centroid.(mm[getindex.(qidx,1)])))
    t2 = time()
    @show t2-t1
    cc = get_num_spikes(vpvrp, vpoc)
    nspikes = round.(Int64, cc)
    nspikes, mpos, mgaze
end

function fit_glm_alt(vpvrp::ViewAndPlaceRepresentationNew, vpoc::ViewAndPlaceOccupancy, unity_gaze_data::UnityRaytraceData)
    # find place bins with minimum amount of occpuancy
    goodbinidx = findall(dropdims(sum(dropdims(sum(vpoc.weight_view,dims=1),dims=1) .> 0.02,dims=2),dims=2).> 5)
    # get the speed per place bin
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, vpoc.placebin_idx,nelements(m_floor))
    cc = get_num_spikes(vpvrp, vpoc)
    nt = size(vpoc.weight_view,3)
    # grab all non-zeros bins 
    qidx1 = vpoc.weight_view[:,goodbinidx,:].>0
    qidx2 = reshape(vv[goodbinidx,:],1,length(goodbinidx), size(vv,2)).>1.0
    qidx = findall(qidx1.&qidx2)
    # replace with goodbinidx 
    qidx2 = [CartesianIndex(ii.I[1], goodbinidx[ii.I[2]], ii.I[3]) for ii in qidx]
    np = length(qidx)
    # FIXME: This is a bit clunky
    nspikes = cc[:,goodbinidx,:][qidx]

    # also grab the correspoding positons and gaze
    # vpoc contains the view and place bin; we need to grab the actual place and gaze from unity_gaze_data
    pos = zeros(2,np)
    gaze = zeros(3,np)
    hd = zeros(np)
    offset = 0
    for tidx in 1:nt 
        _qidx = qidx2[getindex.(qidx2,3) .== nt]

        pidx = getindex.(_qidx, 2)
        vidx = getindex.(_qidx,1)

        trialpos = unity_gaze_data.position[tidx]
        trialhd = unity_gaze_data.head_direction[tidx]
        trialgaze = unity_gaze_data.gaze[tidx]

        poidx = in(pidx).(vpoc.placebin_idx[tidx])
        goidx = in(vidx).(vpoc.viewbin_idx[tidx])

        pgidx = poidx.&goidx
        npg = sum(pgidx)
        pos[:,offset+1:offset+npg] .= trialpos[1:2,pgidx]
        gaze[:,offset+1:offset+npg] .= trialgaze[:,pgidx]
        hd[offset+1:offset+npg] .= trialhd[pgidx]
        offset += npg
    end
    nspikes, pos, gaze, hd
end

"""
Fit GLM models to explain nspikes using pos, head_dir, and gaze as well as combinations thereof
"""
function fit_glm(pos::Matrix{T}, head_dir::Vector{T}, gaze::Matrix{T}, nspikes::Vector{<:Integer}) where T <: Real

    n = size(pos,2)
    # fit pos
    X_pos = permutedims([pos;ones(T, 1, n)])
    q_pos = glm(X_pos, nspikes, Poisson())
    # shuffle to test for significance
    d_pos_shuffled = zeros(1000)
    for i in 1:length(d_pos_shuffled)
        _q = glm(X_pos, shuffle(nspikes), Poisson())
        d_pos_shuffled[i] = deviance(_q)
    end
    is_pos_selective = deviance(q_pos) < percentile(d_pos_shuffled, 1)

    X_gaze = permutedims([gaze;ones(T, 1, n)])
    q_gaze = glm(X_gaze, nspikes, Poisson())
    d_gaze_shuffled = zeros(1000)
    for i in 1:length(d_pos_shuffled)
        _q = glm(X_gaze, shuffle(nspikes), Poisson())
        d_gaze_shuffled[i] = deviance(_q)
    end
    is_gaze_selective = deviance(q_gaze) < percentile(d_gaze_shuffled, 1)

    θ = head_dir
    X_hd = [cos.(θ) sin.(θ) ones(n)]
    q_hd = glm(X_hd, nspikes, Poisson())
    d_hd_shuffled = zeros(1000)
    for i in 1:length(d_pos_shuffled)
        _q = glm(X_hd, shuffle(nspikes), Poisson())
        d_hd_shuffled[i] = deviance(_q)
    end
    is_hd_selective = deviance(q_hd) < percentile(d_hd_shuffled, 1)

    # for combinations,we test the joint fit to the product of the individual fits using AIC
    # pos * gaze
    # joint
    X_pos_gaze = [pos;gaze;ones(T,1,n)]
    q_pos_gaze = glm(permutedims(X_pos_gaze), nspikes, Poisson())
    ll_pos_gaze = loglikelihood(q_pos_gaze)
    aic_pos_gaze = length(q_pos_gaze.pp.beta0) - ll_pos_gaze

    # product
    ll_ind = loglikelihood(q_pos)
    ll_ind += loglikelihood(q_gaze)
    aic_ind = length(q_pos.pp.beta0)+length(q_gaze.pp.beta0) - ll_ind

    # how do we shuffle this in a meaningful way?
    # we can shuffle e.g. view within each place bin for instance
    # compute place ll for data with view shuffled and view for data with place shuffed 
    (deviance_pos = (deviance(q_pos), percentile(d_pos_shuffled, 1)),
     β_pos = q_pos.pp.beta0,
     deviance_gaze = (deviance(q_gaze),percentile(d_gaze_shuffled,1)),
     β_gaze = q_gaze.pp.beta0,
     deviance_hd = (deviance(q_hd), percentile(d_gaze_shuffled, 1)),
     β_hd = q_hd.pp.beta0,
     aic_all = (aic_pos_gaze, aic_ind),
     β_all = q_pos_gaze.pp.beta0)
end

function GLMFit(jocc::JointOccupancy, unity_gaze_data::UnityRaytraceData;redo=false, do_save=true)
    fname = DPHT.filename(GLMFit)
    if !redo && isfile(fname)
        glmfit = load_jld2(GLMFit, fname)
    else
        vpvrp = ViewAndPlaceRepresentationNew()
        nspikes, mpos, mgaze, mhd = fit_glm(vpvrp, jocc, unity_gaze_data)
        res = fit_glm(mpos, mhd, mgaze, nspikes) 
        glmfit = GLMFit(res.β_pos, res.deviance_pos, res.β_gaze, res.deviance_gaze, res.β_hd, res.deviance_hd, res.β_all, res.aic_all)
        if do_save
            save_jld2(glmfit, fname)
        end
    end
    glmfit
end

function do_GLMFit(::Type{T}, celldirs::Vector{String},args...;skip_error=true, kwargs...) where T <: Union{GLMFit, GLMFitH{<:Any}}
    is_gaze_selective = fill(false, length(celldirs))
    is_pos_selective = fill(false, length(celldirs))
    is_hd_selective = fill(false, length(celldirs))
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    @showprogress "Computing GLM fits..." for sessiondir in sessiondirs
        try
            jocc, unity_gaze_data = cd(sessiondir) do
                jocc = JointOccupancy(;kwargs...)
                unity_gaze_data = UnityRaytraceData(;kwargs...)
                jocc, unity_gaze_data
            end
            cidx = findall(allsessiondirs.==sessiondir)
            for (cc,celldir) in zip(cidx,celldirs[cidx])
                glmfit = cd(celldir) do
                    T(args..., jocc, unity_gaze_data;kwargs...)
                end
                if T <: GLMFit
                    is_gaze_selective[cc] = glmfit.deviance_gaze[1] < glmfit.deviance_gaze[2]
                    is_pos_selective[cc] = glmfit.deviance_pos[1] < glmfit.deviance_pos[2]
                    is_hd_selective[cc] = glmfit.deviance_hd[1] < glmfit.deviance_hd[2]
                end
            end
        catch ee
            if skip_error
                @show "Could not proocess $(sessiondir)"
                continue
            else
                rethrow(ee)
           end
       end
    end
    is_pos_selective, is_hd_selective, is_gaze_selective
end

function fit_glm(vpr::ViewRepresentation, vpoc::ViewAndPlaceOccupancy)
    # TODO: Use k-means to cluster 
    mm = vpoc.mm
    # get the Gram matrix
    G = grammatrix(mm)
    # This should be improved to filter out scarce bins etc.
    weight_view = dropdims(sum(vpoc.weight_view,dims=2),dims=2)

    cc = zeros(nelements(mm),Hippocampus.numtrials(vpr))
    for (i,gg) in enumerate(vpr.gaze)
       count_on_manifold!(view(cc,:,i), mm, [Meshes.Point(pg...) for pg in gg])
    end
    widx = findall(weight_view .> 0)
    # do a sub-sample
    tidx = sort(shuffle(1:length(widx))[1:1000])
    X = [centroid(mm[ii.I[1]]) for ii in widx[tidx]]
    # this is a bit clunky; convert from Meshes.Point to matrix
    Xm = zeros(3, length(tidx))
    for i in 1:length(X)
        px = coords(X[i]) 
        Xm[:,i] .= (px.x.val, px.y.val, px.z.val) 
    end
    # TODO: Approximate distance of the data as the distance to the nearest centroid for the beginning
    # and ending point, then shortest distance between the two centroids. This is not going to be completely accurate,
    # but maybe close enough
    Gm = zeros(length(tidx), length(tidx)) 
    for ii in 1:length(tidx)
        for jj in 1:length(tidx)
            Gm[jj,ii] = G[widx[tidx[jj]].I[1], widx[tidx[ii]].I[1]]
        end
    end
    Y = cc[widx[tidx]]
    # FIXME: Gm is singular because we have repeated values
    β = inv(Gm)*Xm'*Y
end

function fit_glm(gaze::AbstractVector{T2},mm::SimpleMesh) where T2 <: AbstractVector{T} where T
    nt = length(gaze)
    # bin counts
    cc = zeros(nelements(mm),nt)
    flat_gaze = T[]
    flat_idx = CartesianIndex{2}[]
    for i in 1:nt
        # assign to element of mm
        if isempty(gaze[i])
            continue
        end
        kidx = mapto(mm, Tuple.(gaze[i]))
        for (j,_kidx) in enumerate(kidx)
            if !isempty(_kidx)
                cc[first(_kidx),i] += 1.0
                push!(flat_gaze, gaze[i][j])
                push!(flat_idx, CartesianIndex(first(_kidx),i))
            end
        end
    end
    cc, flat_gaze, flat_idx
end

function lossfunc2(β::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real},L::AbstractMatrix{<:Real}, α::Real)
    η = X'*β
    λ = exp.(η)
    ll1 = mean(-y.*η + loggamma.(y .+ 1) .+  λ)
    ll2 =α*β'*L*β
    ll1 + ll2
end

function α_penalty!(gg, β::AbstractVector{<:Real}, L::AbstractMatrix{<:Real},α::Real)
    gg .+= α*(L + L')*β
end

function α_penalty!(gg, β::AbstractVector{<:Real}, L::Symmetric{<:Real}, α::Real)
    gg .+= 2*α*L*β
end

function lossfunc2_grad!(gg, β::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real},L::AbstractMatrix{<:Real}, α::Real)
    η = X'*β
    δη = X
    fill!(gg, 0.0)
    ny = length(y)
    @assert size(X,2) == ny
    @assert size(X,1) == length(β)
    for i in eachindex(y) 
        @inbounds yi = y[i]
        @inbounds ηi = η[i]
        eηi = exp(ηi)
        for j in eachindex(β) 
            @inbounds δηji = δη[j,i] 
            @inbounds gg[j] += -yi*δηji + eηi*δηji
        end
    end
    #δll1 = dropdims(mean(-reshape(y,1,length(y)).*δη + δλ,dims=2),dims=2)
    gg ./= ny
    α_penalty!(gg, β, L,α)
end

function lossfunc2_fg!(gg, β::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real},L::AbstractMatrix{<:Real}, α::Real)
end

function lossfunc22(β::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real},L::AbstractMatrix{<:Real}, α::AbstractVector{<:Real})
    η = X'*β
    λ = exp.(η)
    ll1 = mean(-y.*η + loggamma.(y .+ 1) .+  λ)
    ll2 =α.*β'*L*β
    ll1 .+ ll2
end

function lossfunc23(β, X, y,L, α)
    η = X'*β
    λ = exp.(η)
    ll1 = mean(-y.*η + loggamma.(y .+ 1) .+  λ)
    ll2 =α*β'*L*β
    ll1 + ll2
end

function logprob(β, X, y)
    η = X'*β
    λ = exp.(η)
    ll1 = mean(y.*η - loggamma.(y .+ 1) .-  λ)
end

function fit_glm_2(X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real};β0::Union{Nothing,Vector{T}}=nothing,α::T=one(T),show_trace=false,show_progress=false) where T <: Real
    d,n = size(X)
    if β0 === nothing
        β0 = randn(T,d+1)
    elseif size(β0,1) == size(X,1)
        push!(β0, randn(T))
    end
    L2 = zeros(T, size(L,1)+1, size(L,2)+1)
    L2[1:size(L,1), 1:size(L,2)] .= L
    Xq = [X;ones(T, 1, n)]
    lf(β) = lossfunc2(β, Xq, y,L2,α)
    g!(g, β) = lossfunc2_grad!(g, β, Xq, y, L2, α)
    prog = ProgressThresh(1e-8;desc="Minimizing...", enabled=show_progress,showspeed=true)
    function callback(state)
        ProgressMeter.update!(prog, state.g_norm)
        return false
    end
    q = optimize(lf, g!, β0, LBFGS(), Optim.Options(;show_trace=show_trace, show_every=10, callback=callback)) 
end

function get_trainidx(n::Integer, nruns::Integer)
    ntrain = round(Int64, 0.8*n) 
    trainidx = fill(0, ntrain, nruns)
    for r in 1:nruns
        _train_idx = shuffle(1:n)[1:round(Int64, 0.8*n)]
        sort!(_train_idx)
        trainidx[:,r] .= _train_idx
    end
    trainidx
end

function cross_validate(α::AbstractVector{T},X::AbstractMatrix{<:Real}, args...;nruns=10, kwargs...) where T <: Real
    trainidx = get_trainidx(size(X,2),nruns)
    dq = Dict{T,Dict{Symbol,Any}}()
    for _α in α
        β,ll = cross_validate(trainidx, X,args...;α=_α,kwargs...)
        dq[_α] = Dict(:β => β, :ll => ll)
    end
    dq, trainidx
end

"""
Run glm fit on `nruns` separate training and testing sets
"""
function cross_validate(X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real};α::T=one(T),nruns=10,kwargs...) where T <: Real
    trainidx = get_trainidx(size(X,2),nruns)
    cross_validate(trainidx, X, y,L;α=α,kwargs...)
end


function cross_validate(trainidx::Matrix{Int64}, X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real};α::T=one(T),show_progress=false, kwargs...) where T <: Real
    d,n = size(X)
    nruns = size(trainidx,2)
    ll = zeros(nruns)
    β = zeros(d+1,nruns)
    prog = Progress(nruns;desc="Running cross-validation", enabled=show_progress, showspeed=true)
    for r in 1:nruns
        _train_idx = trainidx[:,r]
        sort!(_train_idx)
        test_idx = setdiff(1:n, _train_idx)
        X_train = X[:,_train_idx]
        y_train = y[_train_idx]

        X_test = X[:,test_idx]
        y_test = y[test_idx]

        q = fit_glm_2(X_train, y_train, L;α=α,kwargs...)
        ll[r] = logprob(q.minimizer, [X_test;ones(1,length(test_idx))], y_test)
        β[:,r] .= q.minimizer
        next!(prog)
    end
    β,ll , trainidx
end

function GLMFitH(dims::NTuple{N,Symbol};redo=false, kwargs...) where N
    fname = DPHT.filename(GLMFitH{N}, dims)
    h = process_kwargs(GLMFitH{N};kwargs...)
    if h != 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo && isfile(fname)
        glmfit = load_jld2(GLMFitH{N}, fname)
        if isa(glmfit, JLD2.ReconstructedMutable)
            # missing field
            glmfit = GLMFitH{N}(glmfit.β, glmfit.ll, glmfit.α, glmfit.trainidx, glmfit.nspikes, glmfit.dims, glmfit.qidx, true)
        end
    else 
        jocc, unity_raytrace = cd(DPHT.process_level("session")) do
            jocc = Hippocampus.JointOccupancy(;redo=false)
            ud = Hippocampus.UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
            jocc, ud
        end
        glmfit = GLMFitH(dims, jocc, unity_raytrace;redo=redo, kwargs...)
    end
    glmfit
end

function GLMFitH(dims::NTuple{N,Symbol}, jocc::JointOccupancy, unity_gaze_data::UnityRaytraceData;redo=false, do_save=true,α=10.0.^[-2,-3,-4,-5,-6],nruns=10,show_trace=false,show_progress=false,kwargs...) where N
    fname = DPHT.filename(GLMFitH{N}, dims)
    h = process_kwargs(GLMFitH{N};α=α,nruns=nruns)
    if h != 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    fname_inprogress = replace(fname, ".jld2"=>".jld2.inprogress")
    if !redo && isfile(fname)
        glmfit = load_jld2(GLMFitH{N}, fname)
         if isa(glmfit, JLD2.ReconstructedMutable)
            # missing field
            glmfit = GLMFitH{N}(glmfit.β, glmfit.ll, glmfit.α, glmfit.trainidx, glmfit.nspikes, glmfit.dims, glmfit.qidx, true)
        end
    elseif isfile(fname_inprogress)
        error("$(fname) is currently being computed by another process")
    else
        touch(fname_inprogress)
        # TODO: If we doing joint fit, get the cross-validated alpha from the individua fits first 
        if length(dims) > 1
            use_α = zeros(length(dims))
            validate_α = false
            for (i,d) in enumerate(dims)
                _glmfit = GLMFitH((d,), jocc, unity_gaze_data;α=α,nruns=nruns)
                use_α[i],aidx = get_best_α(_glmfit)
            end
        else
            validate_α = true
            use_α = α
        end
        mm = get_maze_mesh() 
        m_floor = Shadow("xy")(floor_topology3())
        # maybe make this more flexible
        nhd_bins = 24
        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        nspikes, mpos, mgaze, mhd,qidx = fit_glm(vpvrp, jocc, unity_gaze_data)
        # construct X based on dims argument
        d = Int64[] 
        didx = Int64[]
        for dd in dims
            if dd == :g
                push!(d,nelements(mm))
                push!(didx, 1)
            elseif dd == :p
                push!(d, nelements(m_floor))
                push!(didx,2)
            elseif dd == :hd
                push!(d,nhd_bins)
                push!(didx, 3)
            end
        end
        n = length(nspikes)
        X = zeros(sum(d), n)
        # set up laplacian
        L = zeros(sum(d), sum(d))
        offset = 0
        for (j,(nd,dd)) in enumerate(zip(d,dims))
            if dd == :g
                A = adjacencymatrix(mm)
                if validate_α
                    a = 1.0
                else
                    a = use_α[j]
                end
                L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
            elseif dd == :p
                A = adjacencymatrix(m_floor)
                if validate_α
                    a = 1.0
                else
                    a = use_α[j]
                end
                L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
            elseif dd == :hd
                if validate_α
                    a = 1.0
                else
                    a = use_α[j]
                end
                A = get_circular_adjancency(nhd_bins)
                L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
            end
            offset += nd
        end

        for (ii,qq) in enumerate(qidx)
            for (j,_didx) in enumerate(didx)
                offset = sum(d[1:j-1])
                X[offset+qq.I[_didx],ii] = 1.0
            end
        end
        Ls = Symmetric(L)
        if length(dims) == 1
            dq,trainidx = cross_validate(α, X, nspikes, Ls;nruns=nruns,show_trace=show_trace,show_progress=show_progress)
        else
            dq,trainidx = cross_validate(X, nspikes, Ls;nruns=nruns,show_trace=show_trace,show_progress=show_progress)
        end
        β = zeros(size(X,1)+1, nruns, length(α))
        ll = zeros(nruns, length(α))
        for (i,_α) = enumerate(α)
            β[:,:,i] = dq[_α][:β]
            ll[:,i] = dq[_α][:ll]
        end
        glmfit = GLMFitH(β, ll, use_α, trainidx, nspikes, dims, qidx,true)
        if do_save
            save_jld2(glmfit, fname)
        end
        rm(fname_inprogress)
    end
    glmfit
end

"""
Smooth function
"""
function lossfunc4(β::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, α::Matrix{<:Real})
    η = X'*β
    λ = exp.(η)
    ll1 = mean(-y.*η + loggamma.(y .+ 1) .+  λ)
    #ll2 =α*η'*L*η
    ll2 = β'*α*β #smoothing
    ll1 + ll2
end

function get_diff_matrix(α::Vector{T}) where T <: Real
    d = length(α)
    A = zeros(T, d,d+1)
    for i in 1:d
        A[i,i:i+1] = [-α[i],α[i]]
    end
    A
end

function fit_glm_2(X::AbstractMatrix{T}, y::AbstractVector{<:Integer};β0::Union{Nothing,T}=nothing,α::Vector{T}=ones(T,size(X,1)-1)) where T <: Real
    d,n = size(X)
    if β0 === nothing
        β0 = randn(T,d+1)
    end
    A = zeros(d-1,d+1)
    A[:,1:d] = get_diff_matrix(α)
    A2 = A'*A
    Xq = [X;ones(T, 1, n)]
    lf(β) = lossfunc4(β, Xq, y,A2)
    q = optimize(lf, β0, LBFGS();autodiff=:forward) 
end

function fit_glm_3(X::AbstractMatrix{T}, y::AbstractVector{<:Integer},Dm::Matrix{<:T};λ0=1.0, σ0=4.0) where T <: Real
    d,n = size(X)
    Z = λ0.*exp.(-(Dm.^2)./(2*σ0^2))
    zmin = minimum(Z[Z.>0])
    Z[Z.==0] .= zmin
    ll = dropdims(sum(-y.*log.(Z) .+ loggamma.(y .+1 ) .+ Z,dims=1),dims=1)
    ll0,ii = findmin(ll)
    Dm0 = Dm[:,ii]
    # TODO: Refine this further
    # grid search for λ0 and σ 
    # we can actually do optimization here
    function lossfunc3(θ)
        λ,σ = θ
        Z0 = λ.*exp.(-(Dm0.^2)./(2*σ^2))
        ll = sum(-y.*log.(Z0) .+ loggamma.(y .+1 ) .+ Z0)
    end
    q = optimize(lossfunc3, [λ0;σ0])
    (μ_idx=ii, λ=q.minimizer[1], σ=q.minimizer[2], fmin=q.minimum)
end

function lossfunc_lux() end
function fit_glm_lux() end
