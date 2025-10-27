using MAT
using LinearAlgebra
using StatsBase

function get_spatial_maps(celldirs::Vector{String})
    map_adsm = map(celldirs) do celldir
       vpc = MAT.matread(joinpath(celldir, "FiltAll","1px","vmpc.mat"))
       reshape(vpc["vmp"]["data"]["maps_adsm"],40,40)
    end
    cat(map_adsm...,dims=3)
end

function get_view_maps(celldirs::Vector{String})
    map_adsm = map(celldirs) do celldir
       vpc = MAT.matread(joinpath(celldir, "FiltAll","1px","vmsv.mat"))
       vpc["vms"]["data"]["maps_adsm"]
    end
    cat(map_adsm..., dims=3)
end



function get_posterior(f::AbstractArray{T,3},nspikes::Vector{<:Real};τ=one(T)) where T <: Real
    if length(τ) == 1
        _τ = repeat([τ], 1,1, length(nspikes))
    else
        _τ = reshape(τ, 1,1, length(τ))
    end
    #aa = sum([_n.*log.(_f) for (_n,_f) in zip(nspikes, eachslice(f;dims=3))])
    aa = dropdims(sum(permutedims(repeat(nspikes,1,1,1),[3,2,1]).*log.(f),dims=3),dims=3)
    bb = -dropdims(sum(f.*_τ,dims=3),dims=3)
    prb = exp.(aa + bb)
    prb ./= sum(filter(isfinite, prb))
    prb
end

function decode_mean_posterior(prob::Matrix{T}, bins::Tuple{T2,T2}) where T2 <: AbstractVector{<:Real} where T <: Real
    xbins,ybins = bins
    fidx = findall(isfinite, prob)
    fidx[argmax(prob[fidx])]
    μ = zeros(2)
    for i in 1:size(prob,1)
        for j in 1:size(prob,2)
            if isfinite(prob[i,j])
                μ .+= prob[i,j]*[xbins[i], ybins[j]]
            end
        end
    end
    (μ...,)
end

function decode_mean_posterior(prob::Matrix{T}, mm::SimpleMesh) where T <: Real
    _prob = vec(prob)
    fidx = findall(isfinite, _prob)
    μ = zeros(2)
    for ii in fidx
        cp = coords(centroid(mm[ii]))
        μ .+= _prob[ii]*[cp.x.val, cp.y.val]
    end
    μ ./= sum(_prob[fidx])
    (μ...,)
end

function decode_max_posterior(prob::Matrix{T}, bins::Tuple{T2, T2}) where T2 <: AbstractVector{<:Real} where T <: Real
    xbins,ybins = bins
    fidx = findall(isfinite, prob)
    mxi = fidx[argmax(prob[fidx])]
    (xbins[mxi.I[1]], ybins[mxi.I[2]])
end

function decode_max_posterior(prob::Matrix{T}, mm::SimpleMesh) where T <: Real
    d = ndims(mm)
    _prob = vec(prob)
    fidx = findall(isfinite, _prob)
    if isempty(fidx)
        return Tuple(fill(NaN,d))
    end
    mxi = fidx[argmax(_prob[fidx])]
    cp = centroid(mm[mxi])
    Tuple(cp)
end

function decode_patch_posterior(prob::Matrix{T}, bins::Tuple{T2,T2}) where T2 <: AbstractVector{<:Real} where T <: Real
    xbins,ybins = bins
    patches = field_outline(prob;t=1.65)
    qq = [sum(prob[patch]) for patch in patches]
    mxi = argmax(qq)
    μ = zeros(T,2)
    for ii in patches[mxi]
        μ .+= prob[ii].*[xbins[ii.I[1]], ybins[ii.I[2]]]
    end
    μ ./= qq[mxi]
    (μ...,)
end

function decode_patch_posterior(prob::Matrix{T}, mm::SimpleMesh) where T <: Real
    _prob = vec(prob)
    patches = field_outline(_prob,mm;t=1.65)
    qq = [sum(prob[patch]) for patch in patches]
    mxi = argmax(qq)
    d = ndims(mm)
    μ = zeros(T,d)
    for ii in patches[mxi]
        cp = centroid(mm[ii])
        μ .+= _prob[ii].*[Tuple(cp)...,]
    end
    μ ./= qq[mxi]
    (μ...,)
end

function decoded_trajectory(f, nspikes, idx0, idx1, xbins, ybins;τ=0.05)
    nt = size(nspikes,2)
    traj = Vector{Vector{Tuple{Float64, Float64}}}(undef, nt)
    for i in 1:nt
        traj[i] = [decode(get_posterior(f, nspikes[ii,i,:];τ=τ),xbins,ybins)[1] for ii in idx0[i]:min(idx1[i],size(weight,1))]
        #traj[i] = [(xbins[ii.I[1]],ybins[ii.I[2]]) for ii in idx]
    end
    traj
end

function decode_place(sessiondir::String)
    place_selective_cells = open("data/place_selective_cells.txt") do fid
       readlines(fid)
    end

    sessiondirs = DPHT.get_level_path.("session", place_selective_cells)
    cidx = findall(sessiondirs.==sessiondir)
    f = get_spatial_maps(place_selective_cells[cidx])

    # get the points of the place fields
    patches = [Hippocampus.field_outline(f[:,:,i];t=1.65) for i in 1:size(f,3)]
    # create outlines via convex hull
    xbins = range(-12.5f0, stop=12.5f0, length=40);
    ybins = range(-12.5f0, stop=12.5f0, length=40);
    all_hulls = Any[]
    for patch in patches
        _hulls = Any[]
        for p in patch
            points = [Meshes.Point(xbins[ci.I[1]], ybins[ci.I[2]]) for ci in  p]
            chull = convexhull(points)
            push!(_hulls, chull)
        end
        push!(all_hulls, _hulls)
    end
    # find points for with overlaps
    udata = cd(sessiondir) do
        UnityData()
    end
    nt = numtrials(udata)
    spa = map(place_selective_cells[cidx]) do celldir
       cd(celldir) do
        Hippocampus.TrialAlignedSpiketrain()
       end
    end

    weight,bins = Hippocampus.compute_psth(spa, 0.05,8;tmax=20.0)
    hulls = cat(all_hulls..., dims=1)
    timestamps,pos = get_timepoints(udata, hulls)
    idx = CartesianIndex[]
    # find the bin index correpsonding to these time stamps
    allpos = Tuple{Float64,Float64}[]
    for i in 1:nt
        for (pp,tt) in zip(pos[i], timestamps[i])
            kk = searchsortedfirst(bins[1:size(weight,1)], tt)
            if 0 < kk <= size(weight,1)
                push!(idx, CartesianIndex(kk,i))
                push!(allpos, pp)
            end
        end
    end

    XX = zeros(length(idx), size(weight,3))
    for (j,ii) in enumerate(idx)
       XX[j,:] = weight[ii.I[1], ii.I[2],:]
    end

    # grab smoe points to display the error for

    #km_results = kmeans(stack(allpos), 250)
    decoded_pos = Vector{Tuple{Float64, Float64}}(undef, length(allpos))
    for j in axes(XX,1)
        prb = get_posterior(f, XX[j,:])
        μ,fidx = decode(prb, xbins, ybins)
        decoded_pos[j] = μ
        # compute error
    end
    allpos, decoded_pos
end

function plot_fields_with_outline!(ax, f::AbstractArray{T,3}) where T <: Real
    patches = [field_outline(f[:,:,i];t=1.65) for i in 1:size(f,3)]
    xbins = range(-12.5f0, stop=12.5f0, length=40);
    ybins = range(-12.5f0, stop=12.5f0, length=40);
    all_hulls = Any[]
    for patch in patches
        _hulls = Any[]
        for p in patch
            points = [Meshes.Point(xbins[ci.I[1]], ybins[ci.I[2]]) for ci in  p]
            chull = hull(points, JarvisMarch())
            push!(_hulls, chull)
        end
        push!(all_hulls, _hulls)
    end

    colors = to_colormap(:tab20)
    for (i,hulls) in enumerate(all_hulls)
        for _hull in hulls
            b = boundary(_hull)
            if b !== nothing
                viz!(b, color=colors[mod(i-1,10)+1])
            end
        end
    end
end

function plot_fields(patches::Vector{Vector{Vector{CartesianIndex{2}}}})
    colors = to_colormap(:tab10)
    with_theme(plot_theme) do
        fig = Figure()
        ax = Axis(fig[1,1])
        for (ii,patch) in enumerate(patches)
            color = colors[mod(ii-1,10)+1] 
            for p in patch
                plot_field!(ax, p;color=color)
            end
        end
        fig
    end
end

function plot_field!(ax, patch::Vector{CartesianIndex{2}};kwargs...)
    xbins = range(-12.5f0, stop=12.5f0, length=40);
    ybins = range(-12.5f0, stop=12.5f0, length=40);
    points = [Meshes.Point(xbins[ci.I[1]], ybins[ci.I[2]]) for ci in  patch]
    chull = hull(points, GrahamScan())

    b = boundary(chull)
    if b !== nothing
        viz!(ax, b;kwargs...)
    end
end

function plot_field!(ax, patch::Vector{Int64},mm::SimpleMesh;kwargs...)
    hulls = Any[]
    points = [coords(centroid(mm[ci])) for ci in  patch]
    chull = hull(points, JarvisMarch())

    b = boundary(chull)
    if b !== nothing
        viz!(ax, b;kwargs...)
    end
end

function generate_pseudotrial(X::Matrix{T}, assignments::Vector{<:Integer},k::Integer, ntrials::Integer=1) where T <: Real
    Xn = zeros(T,size(X,1),ntrials)
    trialidx = zeros(Int64, size(Xn)...)
    generate_pseudotrial!(Xn, trialidx, X, assignments, k)
    Xn,trialidx
end

function generate_pseudotrial!(Xn::AbstractMatrix{T}, trialidx::AbstractMatrix{Int64}, X::Matrix{T},  assignments::Vector{<:Integer},k::Integer) where T <: Real
    ncells = size(X,1)
    tidx = findall(assignments.==k)
    for j in axes(Xn,2)
        for i in 1:ncells
            ti = rand(tidx)
            Xn[i,j] = X[i,ti]
            trialidx[i,j] = ti
        end
    end
    Xn, trialidx
end

function prepare_decoding_variables(X,Y,twin,tidx=1:size(X,2);tmax=0.2)
    y = zeros(size(Y,1), length(tidx))
    nspikes = zeros(size(X,1), length(tidx))
    tw = zeros(length(tidx))
    bcounts = zeros(length(tidx))
    for i in 1:length(tidx)
       y[:,i],tw[i],didx = Hippocampus.merge_by_time(Y, twin,tidx[i];tmax=tmax)
       nspikes[:,i] = dropdims(sum(X[:,didx],dims=2),dims=2)
       bcounts[i] = length(didx)
    end
    nspikes, y,tw, bcounts
end

function decode_view_simple(X,Y,twin,f,tidx=1:size(X,2))
    maze_points,maze_idx= Hippocampus.map_from_matlab()
    decoded_pos = zeros(3, length(tidx))
    actual_pos = zeros(3, length(tidx))
    for i in 1:length(tidx)
       y,t,didx = Hippocampus.merge_by_time(Y, twin,tidx[i];tmax=0.2)
       nspikes = dropdims(sum(X[:,didx],dims=2),dims=2)
       prb = dropdims(sum(f.*reshape(nspikes,1,1,length(nspikes)),dims=3),dims=3)
       prb ./= sum(filter(isfinite, prb))
       prbf = prb[maze_idx]
       fidx = findall(isfinite.(prbf))
       jk= fidx[argmax(prbf[fidx])]
       decoded_pos[:,i] .= maze_points[jk]
       actual_pos[:,i] .= y
    end
    actual_pos, decoded_pos
end

function decode(m_floor::SimpleMesh, jk)
   cp = coords(centroid(m_floor[jk.I[1]])) 
   (cp.x.val, cp.y.val)
end

function decode(bins::Tuple{AbstractVector{T}, AbstractVector{T}}, jk) where T <: Real
    xbins,ybins = bins
    (xbins[jk.I[1]], ybins[jk.I[2]])
end

function decode_place(X,Y,twin,f,domain;tidx=1:size(X,2), decoder=decode,prog=nothing)
    d = size(Y,1)
    decoded_pos = zeros(d, length(tidx))
    actual_pos = zeros(d, length(tidx))
    if prog === nothing
        prog = Progress(length(tidx),"Decoding...")
    end
    for i in 1:length(tidx)
       y,tw,didx = Hippocampus.merge_by_time(Y, twin,tidx[i];tmax=0.2)
       nspikes = dropdims(sum(X[:,didx],dims=2),dims=2)
       prb = get_posterior(f, nspikes;τ=tw)
       cp = decoder(prb,domain)
       decoded_pos[:,i] .= cp
       actual_pos[:,i] .= y
       next!(prog)
    end
    actual_pos, decoded_pos
end

function decode_place(X::Matrix{<:Real},Y::Matrix{<:Real},twin::AbstractVector{<:Real},f1,f2, domain1, domain2;tidx=1:size(X,2), decoder=decode,prog=nothing)
    d = size(Y,1)
    decoded_pos = zeros(d, length(tidx))
    actual_pos = zeros(d, length(tidx))
    if prog === nothing
        prog = Progress(length(tidx),"Decoding...")
    end
    d1 = embeddim(domain1)
    d2 = embeddim(domain2)
    Y1 = Y[1:d1, :]
    Y2 = Y[d1+1:end,:]
    for i in 1:length(tidx)
       y1,y2,tw,didx = Hippocampus.merge_by_time(Y1, Y2, twin,tidx[i];tmax=0.2)
       nspikes = dropdims(sum(X[:,didx],dims=2),dims=2)
       prb1 = get_posterior(f1, nspikes;τ=tw)
       prb2 = get_posterior(f2, nspikes;τ=tw)
       cp1 = decoder(prb1,domain1)
       cp2 = decoder(prb2,domain2)
       decoded_pos[:,i] .= [cp1...;cp2...]
       actual_pos[1:length(y1),i] .= y1
       actual_pos[length(y1)+1:end,i] .= y2
       next!(prog)
    end
    actual_pos, decoded_pos
end

function compute_place_error_surrogates(X,Y,twin,f,domain, tidx, km_results, decoder=decode;nruns=100)
    prog = Progress(length(tidx)*nruns, "Decoding surrogates...")
    mean_err = zeros(maximum(km_results.assignments), nruns)
    for r in 1:nruns
        qidx = shuffle(1:size(X,2))
        actual_pos, decoded_pos = decode_place(X[:,qidx], Y, twin, f, domain;tidx=tidx, decoder,prog=prog)
        err = sqrt.(dropdims(sum(abs2, decoded_pos .- actual_pos,dims=1),dims=1))
        mean_err[:,r] = vec(Hippocampus.merge_responses(reshape(err, 1, length(err)), km_results.assignments, km_results.counts))
    end
    mean_err
end

function compute_place_error_surrogates(X::Matrix{<:Real},Y::Matrix{<:Real},twin::AbstractVector{<:Real},f1,f2,domain1,domain2, tidx, assignments1, assignments2, decoder=decode;nruns=100)
    prog = Progress(length(tidx)*nruns, "Decoding surrogates...")
    mean_err = fill(NaN, maximum(assignments1), maximum(assignments2),nruns)
    for r in 1:nruns
        qidx = shuffle(1:size(X,2))
        actual_pos, decoded_pos = decode_place(X[:,qidx], Y, twin, f1, f2, domain1,domain2, tidx, decoder;prog=prog)
        err = sqrt.(dropdims(sum(abs2, decoded_pos .- actual_pos,dims=1),dims=1))
        mean_err[:,:,r] = Hippocampus.merge_responses(reshape(err, 1, length(err)), [assignments1,assignments2])
    end
    mean_err
end


function decode_place_simple(X,Y,twin,f,tidx=1:size(X,2))
    m_floor = Shadow("xy")(floor_topology3())
    decoded_pos = zeros(2, length(tidx))
    actual_pos = zeros(2, length(tidx))
    for i in 1:length(tidx)
       y,tw,didx = Hippocampus.merge_by_time(Y, twin,tidx[i];tmax=0.2)
       nspikes = dropdims(sum(X[:,didx],dims=2),dims=2)
       prb = dropdims(sum(f.*reshape(nspikes,1,1,length(nspikes)),dims=3),dims=3)
       fidx = findall(isfinite.(prb))
       jk= fidx[argmax(prb[fidx])]
       cp = coords(centroid(m_floor[jk.I[1]]))

       decoded_pos[:,i] = [cp.x.val, cp.y.val]
       actual_pos[:,i] .= y
    end
    actual_pos, decoded_pos
end

function decode(f::AbstractArray{<:Real,3}, nspikes::Matrix{T},twin::Vector{T}=ones(T, size(nspikes,2))) where T <: Real
    m_floor = Shadow("xy")(floor_topology3())
    decoded_pos = zeros(T, 2, size(nspikes,2))
    @show "here"
    for i in axes(nspikes,2)
       prb = dropdims(sum(f.*reshape(nspikes[:,i],1,1,length(nspikes[:,i])).*twin[i],dims=3),dims=3)
       fidx = findall(isfinite.(prb))
       jk= fidx[argmax(prb[fidx])]
       cp = coords(centroid(m_floor[jk.I[1]]))
       decoded_pos[:,i] = [cp.x.val, cp.y.val]
    end
    decoded_pos
end

function compute_error(point1::Meshes.Point, point2::NTuple{2,T}) where T <: Real
    pc = coords(point1)
    sqrt((point2[1] - pc.x.val)^2 + (point2[2] - pc.y.val)^2)
end 

"""
A summary of the various ways in which we can try and decode space
"""
function plot_spatial_decoding_analysis(;redo=false)
    fname = joinpath(@__DIR__,"..","data","spatial_decoding_analysis.jld2")
    m_floor = Shadow("xy")(Hippocampus.floor_topology3())
    xbins = range(-12.5f0, stop=12.5f0, length=40)
    ybins = range(-12.5f0, stop=12.5f0, length=40)
    if redo || !isfile(fname)

        # get the data
        place_selective_cells = open("data/place_selective_cells.txt") do fid
        readlines(fid)
        end
        spr = map(place_selective_cells) do celldir
            cd(celldir) do
                Hippocampus.SpatialRepresentation(;min_speed=2.0,trial_start=2)
            end
        end

        # combine all spatial responses
        X,Y,twin = Hippocampus.get_population_representation(spr)

        spm = map(place_selective_cells) do celldir
            cd(celldir) do
                Hippocampus.SpatialMapNew(;min_speed=2.0,trial_start=2)
            end
        end

        # smooth all maps with a gaussian filter with 2 bin standard deviation
        spm_smooth = Hippocampus.SmoothedMap.(spm; σ=2,m=5,edge_correct=false)
        f_smooth = cat(Hippocampus.get_rate_map.(spm_smooth)...,dims=3)

        # grab 1000 random data points
        tidx = sort(shuffle(1:size(X,2))[1:1000]);    
        actual_pos, decoded_pos_mean_posterior = Hippocampus.decode_place(X,Y,twin,f_smooth, m_floor, tidx, Hippocampus.decode_mean_posterior)
        err_mean_posterior = sqrt.(dropdims(sum(abs2, decoded_pos_mean_posterior .- actual_pos,dims=1),dims=1))
        actual_pos, decoded_pos_patch_posterior = Hippocampus.decode_place(X,Y,twin,f_smooth, m_floor, tidx, Hippocampus.decode_patch_posterior)
        err_patch_posterior = sqrt.(dropdims(sum(abs2, decoded_pos_patch_posterior .- actual_pos,dims=1),dims=1))
        actual_pos, decoded_pos_max_posterior = Hippocampus.decode_place(X,Y,twin,f_smooth, m_floor, tidx, Hippocampus.decode_max_posterior)
        err_max_posterior = sqrt.(dropdims(sum(abs2, decoded_pos_max_posterior .- actual_pos,dims=1),dims=1))
        JLD2.save(fname, Dict("actual_pos"=>actual_pos, "decoded_pos_mean_posterior"=>decoded_pos_mean_posterior,
                         "decoded_pos_max_posterior"=>decoded_pos_max_posterior,
                         "decoded_pos_patch_posterior"=>decoded_pos_patch_posterior,
                         "err_mean_posterior" => err_mean_posterior,
                         "err_max_posterior" => err_max_posterior,
                         "err_patch_posterior" => err_patch_posterior,
                         "X"=>X,"Y"=>Y, "twin"=>twin,"tidx"=>tidx,"f_smooth"=>f_smooth))
    else
        actual_pos, decoded_pos_max_posterior, decoded_pos_mean_posterior, decoded_pos_patch_posterior = JLD2.load(fname,
                                                                "actual_pos","decoded_pos_max_posterior",
                                                                "decoded_pos_mean_posterior","decoded_pos_patch_posterior")
        err_mean_posterior, err_max_posterior, err_patch_posterior = JLD2.load(fname, "err_mean_posterior",
                                                                               "err_max_posterior",
                                                                               "err_patch_posterior")
    end

    # plot the results
    with_theme(plot_theme) do
        fig = Figure(size=(1024,700))
        lg1 = GridLayout(fig[1,1])
        ax1 = Axis(lg1[1,1])
        viz!(ax1, m_floor;color=:lightgray)
        sc1 = scatter!(ax1, Point2f.(eachcol(decoded_pos_mean_posterior)), color=err_mean_posterior)
        Colorbar(lg1[1,2], sc1, label="Mean posterior error")

        lg2 = GridLayout(fig[1,2])
        ax2 = Axis(lg2[1,1])
        viz!(ax2, m_floor;color=:lightgray)
        sc2 = scatter!(ax2, Point2f.(eachcol(decoded_pos_max_posterior)), color=err_max_posterior)
        Colorbar(lg2[1,2], sc2, label="Max posterior error")

        lg3 = GridLayout(fig[1,3])
        ax3 = Axis(lg3[1,1])
        viz!(ax3, m_floor;color=:lightgray)
        sc3 = scatter!(ax3, Point2f.(eachcol(decoded_pos_patch_posterior)), color=err_patch_posterior)
        Colorbar(lg3[1,2], sc3, label="Patch posterior error")

        lg4 = GridLayout(fig[2,1])
        ax4 = Axis(lg4[1,1])
        viz!(ax4, m_floor;color=:lightgray)
        sc4 = scatter!(ax4, Point2f.(eachcol(actual_pos)), color=:black)

        # plot the field outlines
        lg5 = GridLayout(fig[3,1])
        ax5 = Axis(lg5[1,1])
        viz!(ax5, m_floor;color=:lightgray)
        for i in 1:size(f_smooth,3)
            patches = field_outline(f_smooth[:,:,i],m_floor)
            for patch in patches
                plot_field!(ax, patch,m_floor)
            end
        end
        fig
    end
end

function plot_spatial_decoding_results(km_results, mean_err::Vector{T}, mean_err_sh::Matrix{T};_plot_theme=plot_theme) where T <: Real
    with_theme(_plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_spatial_decoding_results!(lg, km_results, mean_err, mean_err_sh;_plot_theme=_plot_theme)
        fig
    end
end

function plot_spatial_decoding_results!(lg, km_results, mean_err::Vector{T}, mean_err_sh::Matrix{T};_plot_theme=plot_theme) where T <: Real
    err_zscore = (mean_err .- dropdims(mean(mean_err_sh,dims=2),dims=2))./dropdims(std(mean_err_sh, dims=2),dims=2)
    m_floor = Shadow("xy")(Hippocampus.floor_topology3())
    midx = err_zscore .< -2.0
    with_theme(_plot_theme) do
        ax = Axis(lg[1,1],aspect=1)
        viz!(ax, m_floor;color=:lightgray)
        scatter!(ax, Point2f.(eachcol(km_results.centers[:,midx])), color=:red, markersize=12px)
        sc = scatter!(ax, Point2f.(eachcol(km_results.centers)), color=err_zscore, markersize=10px)
        Colorbar(lg[1,2], sc, label="Z-scored error")
    end
end

function plot_view_decoding_results(km_results, mean_err::Vector{T}, mean_err_sh::Matrix{T};_plot_theme=plot_theme) where T <: Real
    with_theme(_plot_theme) do
        fig = Figure(size=(556,538))
        lg = GridLayout(fig[1,1])
        plot_view_decoding_results!(lg, km_results, mean_err, mean_err_sh;_plot_theme=_plot_theme)
        fig
    end
end

function plot_view_decoding_results!(lg, km_results, mean_err::Vector{T}, mean_err_sh::Matrix{T};_plot_theme=plot_theme) where T <: Real
    zscored_err = (mean_err .- dropdims(mean(mean_err_sh,dims=2),dims=2))./dropdims(std(mean_err_sh, dims=2),dims=2)
    mm = get_maze_mesh(;nrefinements=1)
    markersize = fill(7.0, length(zscored_err))
    markersize[zscored_err .< -2.0] .= sqrt(2)*7
    with_theme(_plot_theme) do
        lscene = LScene(lg[1,1],show_axis=false)
        plotmesh!(lscene,mm;segmentcolor=:lightgray, showsegments=true, alpha=0.0, ceiling_offset=10, floor_offset=-20)
        sc = scatter!(lscene, km_results.centers,mm, color=zscored_err, markersize=markersize)
        Colorbar(lg[1,2],sc, label="Z-scored error")
        lscene
    end
end