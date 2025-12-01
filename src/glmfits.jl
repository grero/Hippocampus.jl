using RecurrentNetworkModels
using RecurrentNetworkModels: Lux, Reactant
using RecurrentNetworkModels.Lux: Training
using RecurrentNetworkModels.Optimisers: Adam
using RecurrentNetworkModels:MLDataDevices
using SpecialFunctions
using Random
using Optim
using GLM

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
    @show t2-t1
    cc = get_num_spikes(vpvrp, jocc)
    nspikes = zeros(Int16, length(qidx))
    @show length(intersect(qidx, keys(cc)))
    for (ii,k) in enumerate(qidx)
        if k in keys(cc)
            nspikes[ii] = cc[k]
        end
    end
    nspikes, mpos, mgaze,mhd
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

function GLMFit(celldirs::Vector{String};kwargs...)
    is_gaze_selective = fill(false, length(celldirs))
    is_pos_selective = fill(false, length(celldirs))
    is_hd_selective = fill(false, length(celldirs))
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    @showprogress "Computing GLM fits..." for sessiondir in sessiondirs
        jocc, unity_gaze_data = cd(sessiondir) do
            jocc = JointOccupancy(;kwargs...)
            unity_gaze_data = UnityRaytraceData(;kwargs...)
            jocc, unity_gaze_data
        end
        cidx = findall(allsessiondirs.==sessiondir)
        for (cc,celldir) in zip(cidx,celldirs[cidx])
            glmfit = cd(celldir) do
                GLMFit(jocc, unity_gaze_data)
            end
            is_gaze_selective[cc] = glmfit.deviance_gaze[1] < glmfit.deviance_gaze[2]
            is_pos_selective[cc] = glmfit.deviance_pos[1] < glmfit.deviance_pos[2]
            is_hd_selective[cc] = glmfit.deviance_hd[1] < glmfit.deviance_hd[2]
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
    ll2 =α*η'*L*η
    ll1 + ll2
end

function fit_glm_2(X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real};β0::Union{Nothing,T}=nothing,α::T=T(0.01)) where T <: Real
    d,n = size(X)
    if β0 === nothing
        β0 = randn(T,d)
    end
    lf(β) = lossfunc2(β, X, y,L,α)
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


function lossfunc(model, ps, st, (x,y,L))
    η,st_new = model(x, ps, st)
    ll = sum(RecurrentNetworkModels.poisson_loss.(η,y))
    @show typeof(η) typeof(L)
    ll2 = 0.01f0*η*L*η'
    ll = ll + ll2
    return ll, st_new, (;y_pred=η)
end 

function lossfunc(η, y)
    ll = sum(RecurrentNetworkModels.poisson_loss.(η,y))
    return ll
end

"""
GLM fit using Lux
"""
function fit_glm(X::Matrix{Float32}, y::Vector{<:Integer},L::Matrix{Float32};nepochs=100,learning_rate=Float32(1e-3))

    dev = MLDataDevices.reactant_device()
    cdev = MLDataDevices.cpu_device()

    d,n = size(X)
    model = Lux.Dense(d,1)
    rng = Random.default_rng()
    _ps,_st = Lux.setup(rng, model)
    ps,st = dev((_ps, _st))
    train_state = Training.TrainState(model, ps, st, Adam(learning_rate))

    #mini batch
 

    # validation set
    nval = round(Int64, 0.2*n) 
    validx = shuffle(1:n)[1:nval]
    trainidx = setdiff(1:n, validx)
    xv = X[:,validx]
    yv = reshape(y[validx], 1,nval)
    (xe,ye,Le) = dev.((xv,yv,L[validx,validx]))
    model_compiled = Reactant.@compile model(xe, ps, Lux.testmode(st))

    sort!(validx)
    nt = round(Int64, 0.5*(n-nval))
    ll = lossfunc(model_compiled, ps, st, (xe,ye,Le))
    @show ll[1]
    prog = Progress(nepochs, "Training...")
    #lossfunc = RecurrentNetworkModels.LPoissonLoss()
    for i in 1:nepochs
        tidx = shuffle(trainidx)[1:nt]
        xt = X[:,tidx]
        yt = reshape(y[tidx], 1, length(tidx))
        (xtd,ytd,Ld) = dev.((xt, yt,L[tidx,tidx]))
         (_, loss, _, train_state) = Training.single_train_step!(
                Training.AutoEnzyme(), lossfunc, (xtd, ytd,Ld), train_state
            )
        st_test= Lux.testmode(train_state.states) #
        η_,_ = model_compiled(xe, train_state.parameters, st_test)
        η = cdev(η_)
        lval = lossfunc(η,yv)
        next!(prog;showvalues=[(:lval, lval)])
    end
    model, cdev((train_state.parameters, train_state.states))
end
