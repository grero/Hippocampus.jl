using MDAFiles
using MakieCore
import Makie.SpecApi as S
using StatsBase

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

function extract_spikes(mdata::MountainSortResult, hdata::NeuralData)
    nchs,npts,ntemplates = size(mdata.templates)
    waveforms = Vector{Array{Float64,3}}(undef, ntemplates)
    for i in 1:ntemplates
        template = mdata.templates[:,:,i]
        mx = argmax(abs.(template),dims=2)
        # find the maximum across channels
        imx = argmax(abs.(template[mx]))
        trigger_pt = mx[imx][2]
        idx = findall(mdata.firing[3,:].==i)
        nspikes = length(idx)
        wf = fill(0.0, nchs, npts,nspikes)
        for (j,_idx) in enumerate(idx)
            q = round(Int64,mdata.firing[2,_idx])
            wf[:,:,j] = hdata.data[:,q-trigger_pt+1:q+npts-trigger_pt]
        end
        waveforms[i] = wf
    end
    waveforms
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, x::MountainSortResult)
    # one axis for the templates, another grid layout for auto-correlation
    [PlotSpec(Lines, x.templates[1,:,i]) for i in 1:size(x.templates,3)]
end

function MakieCore.convert_arguments(::Type{<:Lines}, y::Matrix{T}) where T <: Real
    # join each row using NaN
    _y = cat(y,fill(NaN,1,size(y,2)),dims=1)[:]
    _x = repeat([1:(size(y,1)+1);], size(y,2))
    PlotSpec(Lines, _x, _y)
end

function MakieCore.convert_arguments(::Type{<:Plot{plot}}, x::MountainSortResult, y::NeuralData)
    waveforms = extract_spikes(x, y)
    # one axis for each template
    n_spikes_plot = [min(1000, size(wf,3)) for wf in waveforms]
    axes = [S.Axis(plots=[S.Lines(waveforms[i][1,:,1:n_spikes_plot[i]])];xticklabelsvisible=false) for i in 1:length(waveforms)]
    return S.GridLayout(axes)
end