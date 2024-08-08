using MDAFiles
using MakieCore
import Makie.SpecApi as S
using StatsBase
using MultivariateStats
using LinearAlgebra

struct MountainSortResult{T<:Real}
    templates::Array{T,3}
    firing::Matrix{Float64}
end

function discard_templates(mdata::MountainSortResult, cidx)
    func = (!).(in(cidx))
    keep = setdiff(1:size(mdata.templates,3), cidx)
    sort!(keep)
    uidx = fill(0, size(mdata.templates,3))
    uidx[keep] = 1:length(keep)
    to_keep = func.(mdata.firing[3,:])
    new_firing = mdata.firing[:,to_keep]
    # re-index
    new_firing[3,:] = uidx[round.(Int64,new_firing[3,:])]
    new_templates = mdata.templates[:,:,keep]
    MountainSortResult(new_templates, new_firing)

end

function merge_templates(mdata::MountainSortResult, cidx)
    nspikes = countmap(firing[3,:])
    new_template = fill(0.0, size(mdata.templates)[1:2]...)
    new_firing = fill!(similar(mdata.firing), 0.0)
    nn = 0
    for _idx in cidx
       new_template .+= nspikes[_idx]*mdata.templates[:,:,_idx]
       nn += nspikes[_idx]
    end
    new_template ./= nn
    keep = setdiff(1:size(mdata.templates,3), cidx)
    n_cidx = length(keep)+1
    new_templates = cat(mdata.templates[:,:,keep], new_template,dims=3)
    sidx = findall(in(cidx), mda.firing[3,:])
    new_firing[3,sidx] .= n_cidx
    MountainSortResult(new_templates, new_firing)
end


function MountainSortResult()
    templates = MDAFiles.load("templates.mda")
    firings = MDAFiles.load("firings.mda")
    MountainSortResult(templates, firings)
end

function autocorrelogram(spikes::AbstractVector{T},bins::AbstractVector{T};normalize=true) where T <: Real
    nn = length(spikes)
    counts = fill(0.0, length(bins))
    for i in 1:nn-1
        for j in i+1:nn
            Δ = spikes[j]-spikes[i]
            bidx = searchsortedfirst(bins, Δ)
            if bidx <= length(bins)
                counts[bidx] += 1
            else
                break # this works because spiketrains are sorted
            end
        end
    end
    if normalize
        counts ./= sum(counts)
    end
    counts, bins
end

function autocorrelogram(spikes::Vector{T},bins::AbstractVector{T2};kvs...) where T <: AbstractVector{T2} where T2 <: Real
    counts = [fill(0.0, length(bins)) for sp in spikes]
    for (i,sp) in enumerate(spikes)
        counts[i],_ = autocorrelogram(sp,bins;kvs...)
    end
    counts, bins
end

function compute_firing_rate(spikes::Vector{T},windows::AbstractVector{T}) where T <: Real
    λmax = 0.0
    λmean = 0.0
    for i in 1:length(windows)-1
        cc = sum(windows[i] .<= spikes .< windows[i+1])
        cc /= windows[i+1]-windows[1]
        λmax = max(λmax, cc)
        λmean += cc
    end
    λmean /= (length(windows)-1)
    λmax, λmean
end

function compute_firing_rate(spikes::Vector{T},windows::AbstractVector{T2}) where T <: AbstractVector{T2} where T2 <: Real
    nn = length(spikes)
    λmax, λmean = [fill(0.0,nn) for i in 1:2]
    for (i,sp) in enumerate(spikes)
        λmax[i], λmean[i] = compute_firing_rate(sp,windows)
    end
    λmax, λmean
end

function compute_firing_rate(spikes::Vector{T}, ndata::NeuralData;window=1.0) where T <: AbstractVector{T2} where T2 <: Real
    windows = range(0.0, stop=size(ndata.data,2)/ndata.sampling_rate, step=window)
    compute_firing_rate(spikes, windows)
end

function extract_spikes(mdata::MountainSortResult, hdata::NeuralData)
    nchs,npts,ntemplates = size(mdata.templates)
    waveforms = Vector{Array{Float64,3}}(undef, ntemplates)
    spiketrains = Vector{Vector{Float64}}(undef, ntemplates)
    nn = size(hdata.data,2)
    for i in 1:ntemplates
        template = mdata.templates[:,:,i]
        mx = argmax(abs.(template),dims=2)
        # find the maximum across channels
        imx = argmax(abs.(template[mx]))
        trigger_pt = mx[imx][2]
        idx = findall(mdata.firing[3,:].==i)
        _idx = round.(Int64, mdata.firing[2,idx] .+ npts .- trigger_pt) .<= nn
        idx = idx[_idx]
        nspikes = length(idx)
        spiketrain = fill(0.0, nspikes)
        wf = fill(0.0, nchs, npts,nspikes)
        for (j,_idx) in enumerate(idx)
            spiketrain[j] = mdata.firing[2,_idx]/hdata.sampling_rate
            q = round(Int64,mdata.firing[2,_idx])
            wf[:,:,j] = hdata.data[:,q-trigger_pt+1:q+npts-trigger_pt]
        end
        waveforms[i] = wf
        spiketrains[i] = spiketrain
    end
    waveforms,spiketrains
end

"""
Compute PCA on a set of waveforms
"""
function compute_features(::Type{PCA}, waveforms::Matrix{T};nfeatures=3) where T <: Real
    d,n = size(waveforms)
    pca = fit(PCA, waveforms;maxoutdim=nfeatures)
    y = predict(pca, waveforms)
    y
end

function compute_features(::Type{SVD}, waveforms::Matrix{T};nfeatures=3) where T <: Real
    d,n = size(waveforms)
    u,s,v = svd(waveforms)
    permutedims(v[:,1:nfeatures])
end

function compute_features(::Type{TS}, waveforms::Array{T,3};nfeatures=3) where T <: Real where TS
    nch,d,n = size(waveforms)
    features = zeros(T, nch,nfeatures,n)
    for i in axes(features,1)
        features[i,:,:] .= compute_features(TS, waveforms[i,:,:])
    end
    features
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::MountainSortResult)
    # one axis for the templates, another grid layout for auto-correlation
    [PlotSpec(Lines, x.templates[1,:,i]) for i in 1:size(x.templates,3)]
end

Makie.used_attributes(::Matrix{T}) where T <: Real = (:color,)

function MakieCore.convert_arguments(::Type{<:Lines}, y::Matrix{T};color=:black) where T <: Real
    # join each row using NaN
    _y = cat(y,fill(NaN,1,size(y,2)),dims=1)[:]
    _x = repeat([1:(size(y,1)+1);], size(y,2))
    PlotSpec(Lines, _x, _y;color=color)
end

function MakieCore.convert_arguments(::Type{<:Plot{plot}}, x::MountainSortResult, y::NeuralData)
    waveforms,spiketrains = extract_spikes(x, y)
    # look at spike difference up to 500ms with 5 ms step
    bins = range(0.0, step=0.005, stop=0.5)
    counts,_ = autocorrelogram(spiketrains, bins)
    features = compute_features.(SVD, waveforms)
    # firing rate statistics
    λmax, λmean = compute_firing_rate(spiketrains, y)
    nspikes = [size(wf,3) for wf in waveforms]
    if length(features) < 10
        colors=to_colormap(:tab10)
    else
        colors=to_colormap(:tab20)
    end
    # one axis for each template
    n_spikes_plot = min.(1000, nspikes)
    plot_titles = ["$(nspikes[i]) spikes, f_max=$(round(λmax[i], sigdigits=2))Hz, f_mean=$(round(λmean[i],sigdigits=2))Hz" for i in 1:length(nspikes)]
    axes = [S.Axis(plots=[S.Lines(waveforms[i][1,:,1:n_spikes_plot[i]],color=colors[i])];xticklabelsvisible=false,palette=(color=colors,),
                                                                                         ylabel="Unit $i", yticklabelsvisible=false,
                                                                                         topspinevisible=false,
                                                                                         rightspinevisible=false,
                                                                                         xgridvisible=false,
                                                                                         ygridvisible=false,
                                                                                         title=plot_titles[i]) for i in 1:length(waveforms)]

    ax3 = S.Axis3(plots=[S.Scatter(feature[1,1,:], feature[1,2,:],feature[1,3,:];color=color) for (feature,color) in zip(features,colors)];xlabel="PC1",ylabel="PC2",zlabel="PC3",
                                                                                                            xticklabelsvisible=false,
                                                                                                            xticksvisible=false,
                                                                                                            yticklabelsvisible=false,
                                                                                                            yticksvisible=false,
                                                                                                            zticklabelsvisible=false,
                                                                                                            zticksvisible=false,
                                                                                                            zlabeloffset=15,
                                                                                                            xlabeloffset=15,
                                                                                                            ylabeloffset=15,
                                                                                                            )

    ax4 = S.Axis(plots=[S.Lines(bins, counts[i];color=colors[i]) for i in 1:length(counts)];xlabel="Inter-spike-interval [s]")
    return S.GridLayout([S.GridLayout(axes) S.GridLayout(ax3,ax4;rowsizes=[Relative(0.7), Relative(0.3)])];colsizes=[Relative(0.3), Relative(0.7)])
end