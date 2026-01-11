using SpecialFunctions
using Random
using Optim
using GLM
using ReverseDiff
using ADTypes
using HypothesisTests
using SparseArrays
using Dates

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
    testidx::Matrix{Int64}
    nspikes::Vector{Int16}
    dt::Vector{Float64}
    dims::NTuple{N,Symbol}
    qidx::Vector{CartesianIndex{4}}
    converged::Bool
    nrefinements::Vector{Int64}
end

function GLMFitH{N}(β, ll, α, trainidx, nspikes, dims::NTuple{N,Symbol}, qidx) where N
    GLMFitH{N}(β, ll, α, trainidx, nspikes, dims, qidx, true,[3])
end

function GLMFitH(β, ll, α, trainidx, nspikes, dims::NTuple{N,Symbol}, qidx) where N
    GLMFitH{N}(β, ll, α, trainidx, nspikes, dims, qidx, true, [3])
end

"""
    GLMFitH(β::Array{Float64,3}, ll::Matrix{Float64}, α::Vector{Float64}, trainidx::Matrix{Int64}, nspikes::Vector{Int64}, dt, dims::NTuple{N,Symbol}, qidx, converged, nrefinements) where N

Constructor for missing testidx
"""
function GLMFitH(β::Array{Float64,3}, ll::Matrix{Float64}, α::Vector{Float64}, trainidx::Matrix{Int64}, nspikes::Vector{Int16}, dt, dims::NTuple{N,Symbol}, qidx, converged, nrefinements) where N
    n = length(nspikes)
    ntrain, nruns = size(trainidx)
    ntest = n - ntrain
    testidx = fill(0, ntest, nruns)
    for r in 1:nruns
        testidx[:,r] = setdiff(1:n, trainidx[:,r])
    end
    GLMFitH{N}(β, ll, α, trainidx, testidx, nspikes, dt, dims, qidx, converged, nrefinements)
end


DPHT.filename(::Type{GLMFitH{N}}) where N = "glmfith.jld2"

function DPHT.filename(::Type{GLMFitH{N}},dims::NTuple{N,Symbol}) where N
    nn = join(String.(dims),'_')
    fname = "glmfith_$(nn).jld2"
    fname
end

DPHT.level(::Type{GLMFitH{N}}) where N = "cell"

"""
    process_kwargs(::Type{GLMFitH{N}};α=10.0.^[-2,-3,-4,-5,-6], nruns=10, nrefinements=fill(3, N), kwargs...) where N

Return a hash of the supplied keywords.
"""
function process_kwargs(::Type{GLMFitH{N}};α=10.0.^[-2,-3,-4,-5,-6], nruns=10, nrefinements=fill(3, N), trainidx::Union{Matrix{<:Integer}, Nothing}=nothing, testidx::Union{Matrix{<:Integer}, Nothing}=nothing, kwargs...) where N
    h = UInt32(0)
    h = CRC32c.crc32c(string((:α=>α)),h)
    if trainidx !== nothing
        nruns = size(trainidx,2)
    end
    h = CRC32c.crc32c(string((:nruns=>nruns)),h)
    if nrefinements != fill(3,N) 
        h = CRC32c.crc32c(string((:nrefinements=>nrefinements)),h)
    end
    if trainidx !== nothing
        h = CRC32c.crc32c(string((:trainidx=>trainidx)),h)
    end
    if testidx !== nothing
        h = CRC32c.crc32c(string((:testidx=>testidx)),h)
    end
    min_place_obs = get(kwargs, :min_place_obs, 5)
    if min_place_obs != 5
        h = CRC32c.crc32c(string((:min_place_obs=>min_place_obs)),h)
    end
    min_place_dur = get(kwargs, :min_place_dur, 0.02)
    if min_place_dur != 0.02
        h = CRC32c.crc32c(string((:min_place_dur=>min_place_dur)),h)
    end
    ms = get(kwargs, :min_speed, 1) 
    if ms != 1
        h = CRC32c.crc32c(string((:min_speed=>ms)),h)
    end
    h
end

"""
    find_apαpendable(::Type{GLMFitH{N}};α=10.0.^[-2,-3,-4,-5,-6], kwargs...) where N

Return a GLMFitH{N} object where some subset of the specified \alpha s has already been computed. If
no such object can be found, return nothing.
"""
function find_appendable(::Type{GLMFitH{N}}, dims::NTuple{N,Symbol};α=10.0.^[-2,-3,-4,-5,-6], kwargs...) where N
    fname = DPHT.filename(GLMFitH{N}, dims)
    for i in 1:length(α)-1
        h = process_kwargs(GLMFitH{N};α=α[end-i:end],kwargs...)
        hs = string(h, base=16)
        fname_ = replace(fname, ".jld2"=>"_$(hs).jld2")
        if isfile(fname_) 
            glmfit = load_jld2(GLMFitH{N},fname_)
            return glmfit
        end
    end
    return nothing
end




"""
    logprob(glmfit::GLMFitH{N},trainidx::Matrix{Int64}=glmfit.trainidx;in_sample=false) where N

Compute the log-probability of the spiking activity `glmfit.nspikes` recorded using `glmfit.dt` windows
with the fitted coefficients `glmfit.β`.
"""
function logprob(glmfit::GLMFitH{N},trainidx::Matrix{Int64}=glmfit.trainidx;in_sample=false) where N
    if in_sample
        _trialidx = trainidx
    else
        _trialidx = stack(setdiff.([1:length(glmfit.nspikes)], eachcol(trainidx)))
    end
    α,aidx = get_best_α(glmfit)
    β = glmfit.β[:,:,aidx] 
    nb = size(β,1)-1
    if glmfit.dims[1] == :p
        m = 2
    elseif glmfit.dims[1] == :g
        m = 1
    else
        nb = 24
        m = 3
    end
    nspikes = glmfit.nspikes
    X = zeros(nb+1, length(nspikes))
    X[nb+1,:] .= 1.0
    for (k,qq) in enumerate(glmfit.qidx)
        X[qq.I[m],k] = 1.0
    end
    ll = logprob(β, X, nspikes, glmfit.dt,_trialidx)
    ll
end

function compute_edf(dirs::Vector{String}, ::Type{GLMFitH{N}}, args...;skip_error=false, kwargs...) where N
    nf = Dict{String, Vector{Float64}}()
    @showprogress "Computing EDF..." for d in dirs
        cd(d) do
            try
                nf[d] = compute_edf(GLMFitH{N}, args...;kwargs...)
            catch ee
                if !skip_error
                    rethrow(ee)
                end
            end
        end
    end
    nf
end

function compute_edf(::Type{GLMFitH{N}}, dims::NTuple{N,Symbol}, args...;redo=false, do_save=true, append_tag=true, kwargs...) where N
    h = process_kwargs(GLMFitH{N};kwargs...)
    fname = DPHT.filename(GLMFitH{N}, dims)
    fname = replace(fname, ".jld2"=>"_edf.jld2")
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo && isfile(fname)
        nf = JLD2.load(fname, "nf")
    else
        glmfit = GLMFitH(dims, args...;kwargs...)
        nf = compute_edf(glmfit) 
        if do_save
            metadata = Dict{String,Any}() 
            if append_tag
                tag!(metadata, storepatch=true)
            end
            JLD2.save(fname, Dict("meta"=>metadata, "nf"=>nf))
        end
    end
    nf
end

function compute_edf(glmfit)
    X = get_X(glmfit)
    if length(glmfit.dims) > 1
        α = glmfit.α
        aidx = 1
    else
        α,aidx = get_best_α(glmfit)
    end
    β = glmfit.β[:,:,aidx]
    L = get_laplacian(glmfit.dims, (glmfit.nrefinements...,), (α...,))
    compute_edf(β,X, L,glmfit.trainidx)
end

"""
Compute the effective number of degrees of freedom
"""
function compute_edf(β::Matrix{<:Real}, X::Matrix{<:Real}, L::Matrix{<:Real}, trainidx::Matrix{<:Integer})
    L2 = zeros(size(β,1),size(β,1))
    L2[1:end-1,1:end-1] .= L
    # α should already be embedded in L at this point
    nf = zeros(size(β,2))
    @showprogress "Computing edf..." for (i,_trainidx) in enumerate(eachcol(trainidx))
        _X = sparse(X[:,_trainidx])
         W = exp.(Diagonal(vec(β[:,i]'*_X)))
         #W = Diagonal(vec(β[:,i]'*_X))
         H = _X'*inv(_X*W*_X' + L2)*_X*W
         nf[i] = tr(H)
    end
    nf
end

function compute_llrt(::Type{GLMFitH{N}}, args...;kwargs...) where N
    glmfit = GLMFitH(args...;kwargs...)
    nf = compute_edf(GLMFitH{N}, args...;kwargs...)
    compute_llrt(glmfit, nf)
end

function compute_llrt(glmfit_p::GLMFitH{N}, nf::Vector{<:Real};kwargs...) where N
    ll0_in = logprob(glmfit_p.nspikes, glmfit_p.dt, glmfit_p.trainidx;in_sample=true)
    ll1_in = logprob(glmfit_p, glmfit_p.trainidx;in_sample=true)
    llrt = 2*size(glmfit_p.trainidx,1)*(ll1_in .- ll0_in)
    pv = 1 .- cdf.(Chisq.(nf), llrt)
    pv, llrt
end

function compute_llrt(celldirs::Vector{String}, ::Type{GLMFitH{N}}, args...;skip_error=true, kwargs...) where N
    nn = length(celldirs)
    pv = fill(NaN, nn)
    @showprogress "Computing llrt...." for (i,c) in enumerate(celldirs)
        cd(c) do
            try
                _pv, _llrt = compute_llrt(GLMFitH{N},args...;kwargs...)
                nq = sum(_pv .< 0.05)
                pv[i] = 1 - cdf(Binomial(length(_pv), 0.05),nq)
            catch ee
                if !skip_error
                    rethrow(ee)
                end
            end
        end
    end
    pv
end

function plot_glmfit(glmfit::GLMFitH{N},args...;kwargs...) where N
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_glmfit!(lg, glmfit,args...;kwargs...)
        fig
    end
end

function plot_glmfit!(lg, glmfit_joint::GLMFitH{N}, glmfit::Vector{GLMFitH{1}};use_aic=false) where N
    # TODO: Use AIC instead of log-likelihood directly to account for different number of degres of freedom
    nspikes = glmfit_joint.nspikes
    nruns = size(glmfit_joint.ll,1)
    trainidx = glmfit_joint.trainidx
    ntest = length(nspikes) - size(trainidx,1)
    ll = zeros(nruns,length(glmfit))
    dg_joint = size(glmfit_joint.β,1)
    dg_ind = 0
    lgs = [GridLayout(lg[1,i]) for i in 1:length(glmfit)]
    labels = string.(range('A', length=3*length(lgs)+1,step=1))
    for (i,(_lg,_glmfit)) in enumerate(zip(lgs,glmfit))
        aidx = findfirst(_glmfit.α.==glmfit_joint.α[i])
        β = _glmfit.β[:,:,aidx]
        dg_ind += size(β,1)
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
            ll[j,i] += logprob(β[:,j], X, nspikes[testidx])
        end
        plot_glmfit!(_lg, _glmfit, aidx,labels=labels[(i-1)*3+1:i*3];ylabelvisible=i==1)
    end
    #lg2 = GridLayout(lg[2,1])
    Label(lg[2,1,TopLeft()],labels[end])
    ax1 = Axis(lg[2,1])
    ax2 = Axis(lg[2,2])
    if use_aic
        scatter!(ax1, dg_joint .- ntest*glmfit_joint.ll[:,1], dg_ind .- ntest*ll)
    else
        scatter!(ax1, ll[:,1], glmfit_joint.ll[:,1])
        ablines!(ax1, 0.0, 1.0, linestyle=:dot, color=:black)
        scatter!(ax2, ll[:,2], glmfit_joint.ll[:,1])
        ablines!(ax2, 0.0, 1.0, linestyle=:dot, color=:black)
    end
    rowsize!(lg, 1, Relative(0.8))
    ax1.xlabel = "LL $(glmfit[1].dims[1])"
    ax2.xlabel = "LL $(glmfit[2].dims[1])"
    ax1.ylabel = "LL joint"
end

function plot_glmfit!(lg, glmfit::GLMFitH{N},aidx::Union{Int64, Nothing}=nothing;use_aic=false, labels=["A","B","C"], kwargs...) where N
    #show log-likelihoods for each cross-validation fold and the resulting map
    if aidx === nothing
        α,aidx = get_best_α(glmfit)
    end
    if :g in(glmfit.dims)
        pax = LScene 
        mm = get_maze_mesh(;nrefinements=glmfit.nrefinements[1])
        pax_kwargs = (;)
    else
        pax = Axis
        pax_kwargs = (;)
        mm = floor_topology3(;nrefinements=glmfit.nrefinements[1])
    end
    _dof = size(glmfit.β,1)
    ll = glmfit.ll
    nspikes = glmfit.nspikes
    # establish null likelihood
    _dof0 = 1
    ll0 = logprob(glmfit.nspikes, glmfit.dt, glmfit.trainidx, glmfit.testidx)
    xx = repeat(glmfit.α, 1,size(ll,1))[:]
    yy = permutedims(glmfit.ll)[:]
    points = Point2f[]
    for i in 1:size(ll,1)
        for j in 1:size(ll,2)
            push!(points, Point2f(glmfit.α[j],ll[i,j]))
        end
        push!(points, Point2f(NaN))
    end
    lg2 = lg[2,1]
    ax1 = Axis(lg2[1,1],xscale=log10)
    Label(lg2[1,1,TopLeft()],labels[2])
    lines!(ax1, points)
    scatter!(ax1, xx,yy)
    if get(kwargs, :ylabelvisible, true)
        ax1.ylabel = "LL"
    end
    if get(kwargs, :xlabelvisible, true)
        ax1.xlabel = "Smoothing factor"
    end
    vlines!(ax1, glmfit.α[aidx], color=:black, linestyle=:dot)
    ax2 = Axis(lg2[2,1])
    Label(lg2[2,1,TopLeft()],labels[3])
    if use_aic
        scatter!(ax2, _dof .- ntest*ll[:,aidx], _dof0 .- ntest*ll0)
        ax2.xlabel = "AIC"
        
        if get(kwargs, :ylabelvisible, true)
            ax2.ylabel = "AIC null"
        end
    else
        scatter!(ax2, ll[:,aidx], ll0)
        ax2.xlabel = "LL"
        if get(kwargs, :ylabelvisible, true)
            ax2.ylabel = "LL null"
        end
    end
    ablines!(ax2, 0.0, 1.0, linestyle=:dot, color=:black)
    #ax2.xticklabelrotation = -π/7
    axp = pax(lg[1,1];pax_kwargs...)
    Label(lg[1,1,TopLeft()],labels[1])
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

function plot_ll_analysis(glmfit::GLMFitH{N}) where N
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_ll_analysis!(lg, glmfit)
        fig
    end
end

function plot_ll_analysis!(lg, glmfit::GLMFitH{N}) where N
    # get in and out-of sample loglikelihoods for both null and full model
    ll0_in = logprob(glmfit.nspikes, glmfit.dt, glmfit.trainidx;in_sample=true)
    ll0_out = logprob(glmfit.nspikes, glmfit.dt, glmfit.trainidx;in_sample=false)
    ll1_in = logprob(glmfit, glmfit.trainidx;in_sample=true)
    ll1_out = logprob(glmfit, glmfit.trainidx;in_sample=false)
    ax = Axis(lg[1,2])
    ax1 = Axis(lg[1,1])
    ax2 = Axis(lg[2,2])
    ax3 = Axis(lg[2,1])
    scatter!(ax, ll1_out, ll0_out,color=ll1_out.>ll0_out)
    scatter!(ax2, ll1_out, ll1_in)
    scatter!(ax1, ll0_in, ll0_out)
    scatter!(ax3, ll0_in, ll1_in,color=ll1_in.>ll0_in)
    
    linkyaxes!(ax, ax1)
    linkxaxes!(ax, ax2)

    linkyaxes!(ax3, ax2)
    linkxaxes!(ax3, ax1)
    ablines!(ax, 0.0, 1.0, color=:black, linestyle=:dot)
    ablines!(ax1, 0.0, -1.0, color=:red, linestyle=:dot)
    ablines!(ax1, 0.0, 1.0, color=:black, linestyle=:dot)
    ablines!(ax2, 0.0, -1.0, color=:red, linestyle=:dot)
    ablines!(ax1, 0.0, 1.0, color=:black, linestyle=:dot)
    ablines!(ax3, 0.0, 1.0, color=:black, linestyle=:dot)
    ax2.xlabel = "ll out"
    ax2.yticklabelsvisible = false
    ax3.ylabel = "ll in"
    ax3.xlabel = "ll0 in"
    ax.yticklabelsvisible = false
    ax.xticklabelsvisible = false
    ax1.ylabel = "ll0 out"
    ax1.xticklabelsvisible = false
    for _ax in [ax3,ax2]
        _ax.xticklabelrotation = -π/8
    end
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

function logprob(nspikes, w, nruns::Int64)
    trainidx = get_trainidx(length(nspikes),nruns)
    logprob(nspikes, w, trainidx)
end

function logprob(nspikes, w, trainidx::Matrix{Int64};in_sample=false)
    ll = zeros(size(trainidx,2))
    for i in 1:length(ll)
        _trainidx = filter(x->x>0, trainidx[:,i])
        λ = mean(nspikes[_trainidx]./w[_trainidx])
        if !in_sample
            testidx = setdiff(1:length(nspikes), _trainidx)
        else
            testidx = _trainidx
        end
        y = nspikes[testidx]
        ll[i] = mean(y.*log.(λ.*w[testidx]) - loggamma.(y.+1) - λ.*w[testidx])
    end
    ll
end

function logprob(nspikes, w, trainidx::Matrix{Int64}, testidx::Matrix{Int64})
    ll = zeros(size(testidx,2))
    for i in 1:length(ll)
        _trainidx = filter(x->x>0, trainidx[:,i])
        λ = mean(nspikes[_trainidx]./w[_trainidx])
        _testidx = filter(x->x>0, testidx[:,i])
        y = nspikes[_testidx]
        ll[i] = mean(y.*log.(λ.*w[_testidx]) - loggamma.(y.+1) - λ.*w[_testidx])
    end
    ll
end

function logprob(nspikes, w)
    λ = mean(nspikes./w)
    y = nspikes
    ll = mean(y.*log.(λ.*w) - loggamma.(y.+1) - λ.*w)
end

"""
Get the probability that a null model produces the log-likelihoods found in 
cross-validation fold `idx`.
"""
function get_significance_level(glmfit::GLMFitH{N},idx::Integer;α=0.05,use_aic=false, use_signed_rank=false) where N
    # TODO: Why is everything significant?
    ll = glmfit.ll[:,idx]
    nspikes = glmfit.nspikes
    ll0 = zeros(length(ll))
    _dof0 = 1
    _dof = size(glmfit.β,1)
    ntest = length(glmfit.nspikes) - size(glmfit.trainidx,1)
    for i in 1:length(ll0)
        trainidx = glmfit.trainidx[:,i]
        testidx = setdiff(1:length(nspikes), trainidx)
        P = fit(Poisson,nspikes[trainidx])
        ll0[i] = mean(logpdf.(P, nspikes[testidx]))
    end
    if use_aic
        nn = sum(_dof .- ntest*ll .< _dof0 .- ntest*ll0)
    elseif use_signed_rank
        qq =  SignedRankTest(ll, ll0)
        return qq, pvalue(qq, tail=:right)
    else
        nn = sum(ll .> ll0)
    end
    pv = 1 - cdf(Binomial(length(ll),α),nn)
    nn, pv
end

function get_significance_level(glmfit_joint::GLMFitH{N},glmfit::Vector{GLMFitH{1}};α=0.05,use_aic=false, use_signed_rank=false,kwargs...) where N
    ll_joint = glmfit_joint.ll[:,1]
    pv = zeros(length(glmfit))
    ll = zeros(size(ll_joint,1), length(glmfit))
    for (i,_glmfit) in enumerate(glmfit)
        α, aidx = get_best_α(_glmfit)
        ll[:,i] .= logprob(_glmfit, glmfit_joint.trainidx)
        qq = HypothesisTests.SignedRankTest(ll[:,i], _glmfit.ll[:,aidx])
        pv[i]  = pvalue(qq, tail=:left)
    end
    ll,pv
end

function get_significance_level(dims::NTuple{N,Symbol};kwargs...) where N
    glmfit_joint = GLMFitH(dims;kwargs...)
    glmfit = map(dims) do d
        GLMFitH((d,);kwargs...)
    end
    pv = zeros(length(dims))
    for (i,gg) in enumerate(glmfit)
        _,pv[i] = get_significance_level(gg;kwargs...)
    end
    pv_joint = get_significance_level(glmfit_joint, [glmfit...,];kwargs...)
    (pv_joint=pv_joint, pv=pv)
end

function get_significance_level(dims::Symbol;kwargs...) where N
    glmfit = GLMFitH((dims,);kwargs...)
    get_significance_level(glmfit;kwargs...)
end

function get_significance_level(glmfit::GLMFitH{N};α=0.05,use_aic=false, use_signed_rank=true,kwargs...) where N
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
        αb, aidx = get_best_α(glmfit)
        get_significance_level(glmfit, aidx;α=α,use_aic=use_aic,use_signed_rank=use_signed_rank)
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

function get_goodbins(jocc::JointOccupancy;min_place_duration=0.02, min_place_observations=5)
    nt = length(jocc.index)
    nbins = maximum(getindex.(collect(keys(jocc.weight)),2))
    place_weight = zeros(nbins,nt)
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
    goodbinidx = findall(dropdims(sum(place_weight .> min_place_duration,dims=2),dims=2).> min_place_observations)
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

function fit_glm(vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy,unity_gaze_data::UnityRaytraceData;min_place_dur=0.05, min_place_obs=5,min_speed=1.0, min_gaze_dur=0.02, min_gaze_obs=5,kwargs...)
    # TODO: Make sure that vpvrp and jocc use the same trial start reference
    nt = numtrials(unity_gaze_data)
    t1 = time()
    m_floor = Shadow("xy")(floor_topology3())
    mm = get_maze_mesh()
    m_floor = floor_topology3()
    hdbins = range(0.0, stop=2π,length=24)
    place_weight = zeros(nelements(m_floor),nt)
    gaze_weight = zeros(nelements(mm),nt)
    placebin_idx = Vector{Vector{Int64}}(undef, nt)
    gazebin_idx = Vector{Vector{Int64}}(undef, nt)
    for (k,v) in jocc.weight
        vidx,pidx,hidx,tidx = Tuple(k)
        place_weight[pidx,tidx] += v
        gaze_weight[vidx,tidx] += v
    end
    for i in 1:nt
        aidx = jocc.index[i]
        placebin_idx[i] = fill(0, length(aidx))
        gazebin_idx[i] = fill(0, length(aidx))
        for (j,kk) in enumerate(aidx)
            placebin_idx[i][j] = kk.I[2]
            gazebin_idx[i][j] = kk.I[1]
        end
    end
    goodbinidx = findall(dropdims(sum(place_weight .> min_place_dur,dims=2),dims=2).> min_place_obs)
    goodbinidx_g = findall(dropdims(sum(gaze_weight .> min_gaze_dur,dims=2),dims=2).> min_gaze_obs)
    ff = in(goodbinidx)
    ff_g = in(goodbinidx_g)
    # get the speed per place bin
    vv = Hippocampus.compute_speed(unity_gaze_data.position, unity_gaze_data.timestamps, placebin_idx,nelements(m_floor))
    qidx = CartesianIndex{4}[]
    ww = Float64[]
    for (k,v) in jocc.weight
        vidx = k.I[1]
        pidx = k.I[2]
        tidx = k.I[4]
        if ff(pidx) && ff_g(vidx) && (vv[pidx,tidx]  > min_speed)
            push!(qidx,k)
            push!(ww, v)
        end
    end
    # replace with goodbinidx 
    mpos = stack(ustrip.(Meshes.to.(centroid.(m_floor[getindex.(qidx,2)]))))
    mgaze = stack(ustrip.(Meshes.to.(centroid.(mm[getindex.(qidx,1)]))))
    mhd = hdbins[getindex.(qidx,3)] 
    t2 = time()
    cc = get_num_spikes(vpvrp, jocc)
    nspikes = zeros(Int16, length(qidx))
    for (ii,k) in enumerate(qidx)
        if k in keys(cc)
            nspikes[ii] = cc[k]
        end
    end
    # ensure that we are only using non-zero weights
    fidx = ww.>0
    nspikes[fidx], mpos[:,fidx], mgaze[:,fidx],mhd[fidx],qidx[fidx],ww[fidx]
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

function do_GLMFit(::Type{T}, celldirs::Vector{String},dims::NTuple{N,Symbol}, args...;skip_error=true, redo=false, kwargs...) where T <: Union{GLMFit, GLMFitH{<:Any}} where N
    is_gaze_selective = fill(false, length(celldirs))
    is_pos_selective = fill(false, length(celldirs))
    is_hd_selective = fill(false, length(celldirs))
    allsessiondirs = DPHT.get_level_path.("session", celldirs)
    sessiondirs = unique(allsessiondirs)
    nrefinements, kwargs2 = process_refinements(dims;kwargs...)
    @showprogress "Computing GLM fits..." for sessiondir in sessiondirs
        try
            jocc, unity_gaze_data = cd(sessiondir) do
                jocc = JointOccupancy(;nrefinements=nrefinements, kwargs2...)
                unity_gaze_data = UnityRaytraceData(;kwargs...)
                jocc, unity_gaze_data
            end
            cidx = findall(allsessiondirs.==sessiondir)
            for (cc,celldir) in zip(cidx,celldirs[cidx])
                glmfit = cd(celldir) do
                    T(dims,args..., jocc, unity_gaze_data;redo=redo, kwargs...)
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

"""
    fit_glm_all(;kwargs...)

Fit location, gaze, head direction, and combinations of these variables to the cell in the
current directory.

Cross-validation is first performed separately on location and gaze to determine the optimal smoothin
factor for each. The default is to search across 5 values ranging from 10^(-6) to 10^(-2). The
optimal smoothing factor is then used to re-fit using 20 nruns, again for both location and gaze.
We make sure that the same location/gaze bins are used for training of both location and gaze models.
Finally, a model with both location and gaze is trained, again using the optimal smoothign values found above.
"""
function fit_glm_all(; load_only=false, kwargs...)
    # preload these since we need them for all subsequent analysis
    raytrace_fname="unityfile_eyelink_new.csv"
    if load_only
        args = tuple()
    else
        jocc, unity_raytrace = cd(DPHT.process_level("session")) do
            jocc = Hippocampus.JointOccupancy(;redo=false,nrefinements=(p=3,g=2))
            ud = Hippocampus.UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
            jocc, ud
        end
        vpvrp = ViewAndPlaceRepresentationNew(;redo=fname->false, raytrace_fname=raytrace_fname, kwargs...)
        args = (jocc, unity_raytrace, vpvrp)
    end
    # load the base gaze object first
    g = Hippocampus.GLMFitH((:g,),args...; trial_start=2, show_progress=true, nrefinements=[2],load_only=load_only,raytrace_fname=raytrace_fname, kwargs...)
    if g === nothing
        return (;)
    end
    # get the optimal smoothign factor
    α,_ = Hippocampus.get_best_α(g)
    # re-fit for that single α value with more runs
    gg = Hippocampus.GLMFitH((:g,), args...;trial_start=2, show_progress=true, nrefinements=[2],α=[α],nruns=20, raytrace_fname="unityfile_eyelink_new.csv", load_only=load_only,kwargs...)
    if gg === nothing
        return (;)
    end
    # get the base location object
    p = Hippocampus.GLMFitH((:p,),args...;trial_start=2, show_progress=true, nrefinements=[3],load_only=load_only,kwargs...)
    if p === nothing
        return (;glmfit_g=gg)
    end
    # find the optimal smoothing factor
    α,_ = Hippocampus.get_best_α(p)
    # re-fit using the same training idx as for the gaze object above, using the optimal smoothing 
    # factor for location.
    pp = Hippocampus.GLMFitH((:p,),args...;trial_start=2, show_progress=true, nrefinements=[3],α=[α],trainidx=gg.trainidx, testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if pp === nothing
        return (;glmfit_g=gg)
    end

    # head direction
    h = GLMFitH((:hd,),args...; trial_start=2, show_progress=true, nrefinements=[1], raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if h === nothing
        return (glmfit_p=pp, glmfit_g=gg)
    end
    α,_ = Hippocampus.get_best_α(h)
    hh = Hippocampus.GLMFitH((:hd,),args...;trial_start=2, show_progress=true, nrefinements=[1],α=[α],trainidx=gg.trainidx, testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if hh === nothing
        return (glmfit_p=pp, glmfit_g=gg)
    end

    # fit joint pg using the same training idx
    pg = Hippocampus.GLMFitH((:p,:g),args...;trial_start=2, show_progress=true, nrefinements=[3,2], trainidx=gg.trainidx,testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if pg === nothing
        return (glmfit_p=pp, glmfit_g=gg, glmfit_h=hh)
    end
    ph = Hippocampus.GLMFitH((:p,:hd),args...;trial_start=2, show_progress=true, nrefinements=[3,1], trainidx=gg.trainidx,testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if ph === nothing
        return (glmfit_p=pp, glmfit_g=gg, glmfit_h=hh, glmfit_pg=pg)
    end
    gh = Hippocampus.GLMFitH((:g,:hd),args...;trial_start=2, show_progress=true, nrefinements=[2,1], trainidx=gg.trainidx,testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if gh === nothing
        return (glmfit_p=pp, glmfit_g=gg, glmfit_h=hh, glmfit_pg=pg, glmfit_ph=ph)
    end
    pgh = Hippocampus.GLMFitH((:p,:g,:hd),args...;trial_start=2, show_progress=true, nrefinements=[3,2,1], trainidx=gg.trainidx,testidx=gg.testidx, raytrace_fname="unityfile_eyelink_new.csv",load_only=load_only,kwargs...)
    if pgh === nothing
        return (glmfit_p=pp, glmfit_g=gg, glmfit_h=hh, glmfit_pg=pg, glmfit_ph=ph, glmfit_gh=gh)
    end
    (glmfit_g = gg, glmfit_p = pp, glmfit_h=hh, glmfit_pg=pg, glmfit_ph=ph, glmfit_gh=gh, glmfit_pgh=pgh)
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

function lossfunc2(β::AbstractVector{<:Number}, X::AbstractMatrix{<:Number}, y::AbstractVector{<:Number},L::AbstractMatrix{<:Number}, w::AbstractVector{<:Number}, α::Number)
    η = X'*β
    λ = exp.(η)
    ll1 = mean(-y.*η + loggamma.(y .+ 1) .+  λ)
    ll2 =α*β'*L*β
    ll1 + ll2
end

function α_penalty!(gg, β::AbstractVector{<:Number}, L::AbstractMatrix{<:Number},α::Number)
    gg .+= α*(L + L')*β
end

function α_penalty!(gg, β::AbstractVector{<:Number}, L::Symmetric{<:Number}, α::Number)
    gg .+= 2*α*L*β
end

function lossfunc2_grad!(gg, β::AbstractVector{<:Number}, X::AbstractMatrix{<:Number}, y::AbstractVector{<:Number},L::AbstractMatrix{<:Number}, w::AbstractVector{<:Number}, α::Number)
    η = X'*β
    δη = X
    fill!(gg, 0)
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

function logprob(β, X, y,w)
    η = X'*β
    λ = exp.(η)
    ll1 = mean(y.*η - loggamma.(y .+ 1) .-  λ)
end


function logprob(β::Matrix{Float64}, X, y,w, trainidx::Matrix{Int64})
    ll = zeros(size(trainidx,2))
    for i in 1:length(ll)
        _trainidx = filter(x->x>0, trainidx[:,i])
        ll[i] = logprob(β[:,i], X[:,_trainidx], y[_trainidx], w[_trainidx])
    end
    ll
end

"""
Get the dimensions of the GLM model.
"""
function get_dims(dims::NTuple{N,Symbol}, nrefinements::NTuple{N,<:Integer}) where N
    d = Int64[] 
    didx = Int64[]
    for (nf,dd) in zip(nrefinements,dims)
        if dd == :g
            mm = get_maze_mesh(;nrefinements=nf) 
            push!(d,nelements(mm))
            push!(didx, 1)
        elseif dd == :p
            m_floor = Shadow("xy")(floor_topology3(;nrefinements=nf))
            push!(d, nelements(m_floor))
            push!(didx,2)
        elseif dd == :hd
            push!(d,nhd_bins)
            push!(didx, 3)
        end
    end
    d,didx
end

function fit_glm_2(X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real},w::AbstractVector{<:Real};β0::Union{Nothing,Vector{T}}=nothing,α::T=one(T),show_trace=false,show_progress=false,kwargs...) where T <: Real
    d,n = size(X)
    if β0 === nothing
        β0 = randn(T,d+1)
    end
    if size(β0,1) > size(X,1)
        L2 = zeros(T, size(L,1)+1, size(L,2)+1)
        L2[1:size(L,1), 1:size(L,2)] .= L
        Xq = [X;ones(T, 1, n)]
    else
        L2 = L
        Xq = X
    end
    lf(β) = lossfunc2(β, Xq, y,L2,w, α)
    g!(g, β) = lossfunc2_grad!(g, β, Xq, y, L2, w, α)
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

function get_testidx(n::Integer, trainidx::Matrix{Int64})
    ntrain,nruns = size(trainidx)
    ntest = n - ntrain
    testidx = fill(0, ntest, nruns)
    for r in 1:nruns
        testidx[:,r] = setdiff(1:n, trainidx[:,r])
    end
    testidx
end

function get_train_test_idx(n::Integer, nruns::Integer, nchunks::Integer;use_contiguous_chunks=true)
    # divide into chunks first
    if use_contiguous_chunks
        # split into (roughly) equal chunks
        chunksize = round(Int64, floor(n/nchunks))
        cidx = fill(0,n)
        for i in 1:nchunks
            cidx[(i-1)*chunksize+1:i*chunksize] .= i
        end
        cidx[cidx.==0] .= nchunks 
    else
        cidx = rand(1:chunksize, n)
    end
    trainidx = fill(0, n, nruns)
    testidx = fill(0, n, nruns)
    offset_train = 0
    offset_test = 0
    for c in 1:chunksize
        tidx = findall(cidx.==c)
        m = length(tidx)
        ntrain = round(Int64, 0.8*m)
        ntest = m - ntrain
        for r in 1:nruns
            _trainidx = shuffle(tidx)[1:ntrain]
            sort!(_trainidx)
            trainidx[offset_train+1:offset_train+ntrain,r] .= _trainidx
            _testidx = setdiff(tidx, _trainidx)
            testidx[offset_test+1:offset_test+ntest,r] .= _testidx
        end
        offset_train += ntrain
        offset_test += ntest
    end
    trainidx, testidx
end

function cross_validate(α::AbstractVector{T},X::AbstractMatrix{<:Real}, y,L,w;nruns=10, kwargs...) where T <: Real
    trainidx,testidx = get_train_test_idx(size(X,2),nruns,5)
    cross_validate(α,trainidx, testidx, X, y, L, w;kwargs...)
end

function cross_validate(α::AbstractVector{T},trainidx::Matrix{<:Integer}, X::AbstractMatrix{<:Number}, y,L,w;kwargs...) where T <: Real
    testidx = get_testidx(size(X,2), trainidx)
    cross_validate(α, trainidx, testidx, X, y, L, w;kwargs...)
end

function cross_validate(α::AbstractVector{T},trainidx::Matrix{<:Integer}, testidx::Matrix{<:Integer}, X::AbstractMatrix{<:Number}, y,L,w;kwargs...) where T <: Real
    dq = Dict{T,Dict{Symbol,Any}}()
    # TODO: We need to compile to reactant here
    for _α in α
        β,ll,_,cv = cross_validate(trainidx, testidx, X,y,L,w;α=_α,kwargs...)
        dq[_α] = Dict(:β => β, :ll => ll, :converged => cv)
    end
    dq, trainidx
end

"""
Run glm fit on `nruns` separate training and testing sets
"""
function cross_validate(X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real},ww::AbstractVector{<:Real};α::T=one(T),nruns=10,kwargs...) where T <: Real
    trainidx,testidx = get_train_test_idx(size(X,2),nruns,5)
    cross_validate(trainidx, X, y,L,ww;α=α,kwargs...)
end


function cross_validate(trainidx::Matrix{Int64}, X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real},ww::AbstractVector{<:Real};kwargs...) where T <: Real
    testidx = get_testidx(size(X,2), trainidx)
    cross_validate(trainidx, testidx, X, y, L, ww;kwargs...)
end

function cross_validate(trainidx::Matrix{Int64}, testidx::Matrix{Int64}, X::AbstractMatrix{T}, y::AbstractVector{<:Integer}, L::AbstractMatrix{<:Real},ww::AbstractVector{<:Real};α::T=one(T),show_progress=false, do_confound=false, kwargs...) where T <: Real
    d,n = size(X)
    nruns = size(trainidx,2)
    ll = zeros(nruns)
    β = zeros(d+1,nruns)
    has_converged = fill(false, nruns)
    prog = Progress(nruns;desc="Running cross-validation", enabled=show_progress, showspeed=true)
    for r in 1:nruns
        _train_idx = filter(x->x>0, trainidx[:,r])
        sort!(_train_idx)
        _test_idx = filter(x->x>0, testidx[:,r])
        X_train = X[:,_train_idx]
        if do_confound
            _train_idx = filter(x->x>0, trainidx[:,rand(setdiff(1:nruns, r))])
            y_train = y[_train_idx]
            w_train = ww[_train_idx]
        else
            y_train = y[_train_idx]
            w_train = ww[_train_idx]
        end

        X_test = X[:,_test_idx]
        y_test = y[_test_idx]
        w_test = ww[_test_idx]

        q = fit_glm_2(X_train, y_train, L, w_train;α=α,kwargs...)
        ll[r] = logprob(q.minimizer, [X_test;ones(1,length(_test_idx))], y_test, w_test)
        β[:,r] .= q.minimizer
        has_converged[r] = Optim.converged(q)
        next!(prog)
    end
    β,ll , trainidx, has_converged
end

function process_refinements(dims::NTuple{N,Symbol};kwargs...) where N
    _nrefinements = get(kwargs, :nrefinements, [3])
    kwargs2 = filter(k->k[1]!=:nrefinements, kwargs)
    nrf = [3,2]
    nk = (:p,:g)
    for (p,k) in zip(nk,nrf)
        if p in dims
            ii = findfirst(dims.==p)
            jj = findfirst(nk.==p)
            nrf[jj] = _nrefinements[ii]
        end
    end
    nrefinements = NamedTuple{nk}(nrf) 
    nrefinements, kwargs2
end

function isolderthan(fname, tt=now(UTC))
    st = stat(fname)
    if st.mtime < datetime2unix(tt)
        return true
    end
    return false
end

function GLMFitH(dims::NTuple{N,Symbol};redo::Function=fname->false, do_append=false, kwargs...) where N
    fname = DPHT.filename(GLMFitH{N}, dims)
    h = process_kwargs(GLMFitH{N};kwargs...)
    if h != 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = false
    if !redo(fname) && isfile(fname)
        glmfit = load_jld2(GLMFitH{N}, fname)
        if isa(glmfit, JLD2.ReconstructedMutable)
            # missing field
            # this is super hacky
            if !(:dt in typeof(glmfit).parameters[2])
                do_compute = true
            else
                # we should be able to reconstruct here
                glmfit = GLMFitH(glmfit.β, glmfit.ll, glmfit.α, glmfit.trainidx, glmfit.nspikes, glmfit.dt, glmfit.dims, glmfit.qidx, true,glmfit.nrefinements)
            end
        end
    elseif get(kwargs, :load_only, false)
        # no file exist but we requested load only
        return nothing
    else
        do_compute = true
    end
    if do_compute
        if do_append
            glmfit_a = find_appendable(GLMFitH{N},dims;kwargs...)
        else
            glmfit_a = nothing
        end
        nrefinements,kwargs2 = process_refinements(dims;kwargs...) 
        jocc, unity_raytrace = cd(DPHT.process_level("session")) do
            jocc = Hippocampus.JointOccupancy(;redo=false,nrefinements=nrefinements,kwargs2...)
            ud = Hippocampus.UnityRaytraceData(raytrace_fname="unityfile_eyelink_new.csv";redo=false)
            jocc, ud
        end
        glmfit = GLMFitH(dims, jocc, unity_raytrace;redo=redo, appendto=glmfit_a, kwargs...)
    end
    glmfit
end

function get_laplacian(dims::NTuple{N,Symbol}, nrefinements::NTuple{N,<:Integer}, α::NTuple{N,<:Real}) where N
    d,didx = get_dims(dims, nrefinements)
    L = zeros(sum(d), sum(d))
    offset = 0
    for (j,(nd,dd,nf,a)) in enumerate(zip(d,dims,nrefinements,α))
        if dd == :g
            mm = get_maze_mesh(;nrefinements=nf) 
            A = adjacencymatrix(mm)
            L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
        elseif dd == :p
            m_floor = Shadow("xy")(floor_topology3(;nrefinements=nf))
            A = adjacencymatrix(m_floor)
            L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
        elseif dd == :hd
            A = get_circular_adjancency(nhd_bins)
            L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
        end
        offset += nd
    end
    L
end

function get_X(glmfit::GLMFitH{N}) where N
    get_X(glmfit.qidx, glmfit.dims, (glmfit.nrefinements...,))
end

function get_X(qidx, dims, nrefinements)
    n = length(qidx)
    d,didx = get_dims(dims, nrefinements)
    X = zeros(sum(d)+1, n)
    X[end,:] .= 1.0
    for (ii,qq) in enumerate(qidx)
        for (j,_didx) in enumerate(didx)
            offset = sum(d[1:j-1])
            X[offset+qq.I[_didx],ii] = 1.0
        end
    end
    X
end

function GLMFitH(dims::NTuple{N,Symbol}, jocc::JointOccupancy, unity_gaze_data::UnityRaytraceData,vpvrp::Union{ViewAndPlaceRepresentationNew,Nothing}=nothing;appendto::Union{Nothing, GLMFitH{N}}=nothing, redo::Function=fname->false, do_save=true,load_only=false, α=10.0.^[-2,-3,-4,-5,-6],nruns=10,show_trace=false,show_progress=false,nrefinements=fill(3,length(dims)), trainidx::Union{Matrix{<:Integer},Nothing}=nothing, testidx::Union{Matrix{<:Integer}, Nothing}=nothing, kwargs...) where N
    fname = DPHT.filename(GLMFitH{N}, dims)
    h = process_kwargs(GLMFitH{N};α=α,nruns=nruns, nrefinements=nrefinements,trainidx=trainidx,testidx=testidx)
    if h != 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    fname_inprogress = replace(fname, ".jld2"=>".jld2.inprogress")
    do_compute = false
    if !redo(fname) && isfile(fname)
        glmfit = load_jld2(GLMFitH{N}, fname)
         if isa(glmfit, JLD2.ReconstructedMutable)
            # missing field
            if !(:dt in fieldnames(typeof(glmfit)))
                # redo
                do_compute = true
            else
                glmfit = GLMFitH{N}(glmfit.β, glmfit.ll, glmfit.α, glmfit.trainidx, glmfit.nspikes, glmfit.dt, glmfit.dims, glmfit.qidx, true, [3])
            end
        end
    elseif load_only
        return nothing
    elseif isfile(fname_inprogress) && do_save
        error("$(fname) is currently being computed by another process")
    else
        do_compute = true
    end
    if do_compute
        touch(fname_inprogress)
        if vpvrp === nothing
            vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        end
        nspikes, mpos, mgaze, mhd,qidx,ww = fit_glm(vpvrp, jocc, unity_gaze_data)
        if trainidx === nothing
            trainidx,testidx = get_train_test_idx(length(nspikes),nruns,5)
        else
            if testidx === nothing
                testidx = get_testidx(length(nspikes), trainidx)
            end
        end
        nruns = size(trainidx, 2)
        # TODO: If we are doing joint fit, get the cross-validated alpha from the individua fits first 
        if length(dims) > 1
            use_α = zeros(length(dims))
            validate_α = false
            for (i,(d,rf)) in enumerate(zip(dims,nrefinements))
                # we are only doign this to get the \alpha  values
                _glmfit = GLMFitH((d,), jocc, unity_gaze_data;α=α,nrefinements=[rf],kwargs...)
                use_α[i],aidx = get_best_α(_glmfit)
            end
        else
            validate_α = true
            use_α = α
        end
        # maybe make this more flexible
        nhd_bins = 24
        
        # construct X based on dims argument
        d = Int64[] 
        didx = Int64[]
        for (nf,dd) in zip(nrefinements,dims)
            if dd == :g
                mm = get_maze_mesh(;nrefinements=nf) 
                push!(d,nelements(mm))
                push!(didx, 1)
            elseif dd == :p
                m_floor = Shadow("xy")(floor_topology3(;nrefinements=nf))
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
        for (j,(nd,dd,nf)) in enumerate(zip(d,dims,nrefinements))
            if dd == :g
                mm = get_maze_mesh(;nrefinements=nf) 
                A = adjacencymatrix(mm)
                if validate_α
                    a = 1.0
                else
                    a = use_α[j]
                end
                L[offset+1:offset+nd, offset+1:offset+nd] = a*(diagm(dropdims(sum(A,dims=2),dims=2)) - A)
            elseif dd == :p
                m_floor = Shadow("xy")(floor_topology3(;nrefinements=nf))
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
        
        α0 = use_α
        β0 = zeros(size(X,1)+1, nruns, length(use_α)) 
        ll0 = zeros(nruns,length(use_α))
        idx1 = 1:length(use_α)
        if appendto !== nothing
            # check if there is object we can append to
            # use the same training idx
            trainidx .= appendto.trainidx
            testidx .= appendto.testidx
            α_ = appendto.α
            α0 = sort(unique([α_;use_α]),rev=true)
            idx0 = [findfirst(α0.==a0) for a0 in α_]
            idx1 = setdiff(1:length(α0), idx0)
            ll0 = zeros(nruns, length(α0)) 
            ll0[:,idx0] .= appendto.ll
            β0 = zeros(size(X,1)+1, nruns, length(α0)) 
            β0[:,:,idx0] = appendto.β
        end
        if length(dims) == 1
            dq,_ = cross_validate(α0[idx1], trainidx, testidx, X, nspikes, Ls,ww;show_trace=show_trace,show_progress=show_progress,kwargs...)
            for (i,_α) = enumerate(α0[idx1])
                β0[:,:,idx1[i]] = dq[_α][:β]
                ll0[:,idx1[i]] = dq[_α][:ll]
            end
        else
            β,ll,_ = cross_validate(trainidx, testidx, X, nspikes, Ls,ww;show_trace=show_trace,show_progress=show_progress,kwargs...)
            β = reshape(β, size(β)...,1)
            β0[:,:,idx1] .= β
            ll = reshape(ll, size(ll)...,1)
            ll0[:,idx1] .= ll
        end
        glmfit = GLMFitH(β0, ll0, α0, trainidx, testidx, nspikes, ww, dims, qidx,true,nrefinements)
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
