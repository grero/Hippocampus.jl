# functions to compute conjunctions

struct FieldConjunctions{T<:AbstractResponseFields}
    fields::T
    λ_infield::Vector{Vector{Float64}}
    λ_outfield::Vector{Vector{Float64}}
    sic_infield::Vector{Float64}
    sic_outfield::Vector{Float64}
    sic_sub::Vector{Vector{Float64}}
end

function Hippocampus.issignificant(fj::FieldConjunctions;pv_threshold=0.05)
    res = fill(false, length(fj.λ_infield))
    for ii in 1:length(res)
        res[ii] = fj.sic_infield[ii] > percentile(fj.sic_sub[ii], 100*(1-pv_threshold)) 
    end
    res
end

get_mesh(jf::FieldConjunctions{SpatialResponseFields}) = get_maze_mesh(;nrefinements=jf.fields.args[:nrefinements].g) 
get_mesh(jf::FieldConjunctions{GazeResponseFields}) = Shadow("xy")(floor_topology3(;nrefinements=jf.fields.args[:nrefinements].p))

"""
Condition on either view or place fields and compare firing rates to the unconditioned

Compare the firing rates for e.g. view when the animal was in a particular place field vs 
    when it was elsewhere; if there is indeed a conjunction, the conditioned firing rate should be higher.
"""
function conjunctions(jm::JointMap, fields::SpatialResponseFields)
    m_floor = get_mesh(SpatialResponseFields, fields.args[:nrefinements])
    mm = get_maze_mesh(;nrefinements=fields.args[:nrefinements].g)
    # should we loop over fields?
    # treat each field separately?
    clusters = Hippocampus.merge_fields(m_floor, fields.binidx)
    #covered = Set{Int64}()
    #for field in fields
    #    bidx = findall(Meshes.intersects.(m_floor, field))
    #    for b in bidx
    #        push!(covered, b)
    #    end
    #end
    all_not_covered = setdiff(1:nelements(m_floor), fields.binidx)
    res = Dict()
    for (ii,covered) in enumerate(clusters)
        not_covered = setdiff(1:nelements(m_floor), covered)
        not_covered_sub = [shuffle(all_not_covered)[1:length(covered)] for _ in 1:1000]
        # compare firing rates within vis outside the field
        x_covered = zeros(nelements(mm))
        w_covered = zeros(nelements(mm))
        x_not_covered = zeros(nelements(mm))
        w_not_covered = zeros(nelements(mm))

        x_not_covered_sub = zeros(nelements(mm),length(not_covered_sub))
        w_not_covered_sub = zeros(nelements(mm), length(not_covered_sub))
        n_covered = 0
        for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
            pidx = getindex(qidx,2)
            vidx = getindex(qidx,1)
            if pidx in covered
                #push!(Z_covered, λ)
                n_covered += 1
                x_covered[vidx] += w
                w_covered[vidx] += occ
            elseif pidx in not_covered
                #push!(Z_not_covered, λ)
                x_not_covered[vidx] += w
                w_not_covered[vidx] += occ
                for (jj,nc) in enumerate(not_covered_sub)
                    if pidx in nc
                        x_not_covered_sub[vidx,jj] += w
                        w_not_covered_sub[vidx,jj] += occ
                    end

                end
            end
        end
        x_covered = laplace_smoothing(x_covered, mm, 0.01;niter=100)
        w_covered = laplace_smoothing(w_covered, mm, 0.01;niter=100)
        λ_covered  = x_covered./w_covered
        sic_covered = compute_skaggs_sic(λ_covered[:], w_covered[:])

        x_not_covered = laplace_smoothing(x_not_covered, mm, 0.01;niter=100)
        w_not_covered = laplace_smoothing(w_not_covered, mm, 0.01;niter=100)
        λ_not_covered = x_not_covered./w_not_covered
        sic_not_covered = compute_skaggs_sic(λ_not_covered[:], w_not_covered[:])
        sic_sub = zeros(size(x_not_covered_sub,2))
        for ii in 1:size(x_not_covered_sub,2)
            sic_sub[ii] = compute_skaggs_sic(x_not_covered_sub[:,ii]./w_not_covered_sub[:,ii], w_not_covered_sub[:,ii])
        end
        # create a null distribution for the difference in SIC between in-field and out of field.
        res[ii] = Dict("λ_infield" =>λ_covered, "λ_outfield" =>λ_not_covered, "sic_infield"=>sic_covered,
                       "sic_outfield"=>sic_not_covered, "sic_sub"=>sic_sub)
    end
    res
end

abstract type AbstractFieldConjunctions end

struct PlaceViewConjunction <: AbstractFieldConjunctions
    spatial_fields::SpatialResponseFields
    view_fields::GazeResponseFields
    λ_covered::Matrix{Float64}
    λ_sub::Array{Float64,3}
    λ_infield::Matrix{Float64}
    λ_outfield::Vector{Float64}
end
struct ViewPlaceConjunction <: AbstractFieldConjunctions
    view_fields::GazeResponseFields
    spatial_fields::SpatialResponseFields
    λ_covered::Matrix{Float64}
    λ_sub::Array{Float64,3}
    λ_infield::Matrix{Float64}
    λ_outfield::Vector{Float64}
end

function Hippocampus.issignificant(pvc::AbstractFieldConjunctions;pv_threshold=0.05)
    pv = fill(false, size(pvc.λ_covered)[1:2]...)
    for ii in CartesianIndices(pv) 
        vv = filter(isfinite, pvc.λ_sub[:,ii])
        if !isempty(vv)
            qt = percentile(vv, 100*(1-pv_threshold))
            pv[ii] = pvc.λ_covered[ii] > qt
        end
    end
    pv
end

function conjunctions(jm::JointMap, spatial_fields::SpatialResponseFields, view_fields::GazeResponseFields)
    m_floor = get_mesh(SpatialResponseFields, spatial_fields.args[:nrefinements])
    mm = get_maze_mesh(;nrefinements=spatial_fields.args[:nrefinements].g)

    spatial_clusters = Hippocampus.merge_fields(spatial_fields)
    view_clusters = Hippocampus.merge_fields(view_fields)

    all_not_covered = setdiff(1:nelements(m_floor), spatial_fields.binidx)
    res = Dict()
    λ_covered = zeros(length(view_clusters), length(spatial_clusters))
    λ_sub = zeros(1000, length(view_clusters), length(spatial_clusters))
    for (ii,covered) in enumerate(spatial_clusters)
        not_covered = setdiff(1:nelements(m_floor), covered)
        not_covered_sub = [shuffle(all_not_covered)[1:length(covered)] for _ in 1:1000]
        # compare firing rates within vis outside the field
        x_covered = zeros(nelements(mm))
        w_covered = zeros(nelements(mm))
        x_not_covered = zeros(nelements(mm))
        w_not_covered = zeros(nelements(mm))

        x_not_covered_sub = zeros(nelements(mm),length(not_covered_sub))
        w_not_covered_sub = zeros(nelements(mm), length(not_covered_sub))
        n_covered = 0
        for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
            pidx = getindex(qidx,2)
            vidx = getindex(qidx,1)
            if pidx in covered
                #push!(Z_covered, λ)
                n_covered += 1
                x_covered[vidx] += w
                w_covered[vidx] += occ
            elseif pidx in not_covered
                #push!(Z_not_covered, λ)
                x_not_covered[vidx] += w
                w_not_covered[vidx] += occ
                for (jj,nc) in enumerate(not_covered_sub)
                    if pidx in nc
                        x_not_covered_sub[vidx,jj] += w
                        w_not_covered_sub[vidx,jj] += occ
                    end

                end
            end
        end
        # now we have firing rates for each of view bins conditioned on a particular place field
        # aggregate within each of the view fields
        λ_covered[:,ii]  = [mean(filter(isfinite, x_covered[view_fields.binidx[c]]./w_covered[view_fields.binidx[c]])) for c in view_clusters]
        for (jj,cluster) in enumerate(view_clusters)
            qidx = view_fields.binidx[cluster]
            for kk in 1:1000
                ll = x_not_covered_sub[qidx,kk]./w_not_covered_sub[qidx,kk]
                λ_sub[kk,jj,ii] = mean(filter(isfinite, ll))
            end
        end
    end
    PlaceViewConjunction(spatial_fields, view_fields, λ_covered,λ_sub)
end

function process_kwargs(::Type{PlaceViewConjunction},h::UInt32=zero(UInt32);kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    h = process_kwargs(GazeResponseFields,h;kwargs...)
    h
end

function PlaceViewConjunction(;redo=fname->false, do_save=true, kwargs...)
    fname = "place_view_conjunction.jld2"
    h = process_kwargs(PlaceViewConjunction;kwargs...)
    if h > 0
        hs = string(h,base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        X = load_jld2(PlaceViewConjunction,fname)
    else
        jm = JointMap(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields,10_000;kwargs...)
        rf_gaze = get_response_fields(GazeResponseFields,10_000;kwargs...)
        X = conjunctions(jm, rf_spatial, rf_gaze)
        if do_save
            save_jld2(X,fname)
        end
    end
    X
end

function conjunctions(::Type{SpatialResponseFields}, celldir::String;kwargs...)
    jm,rf = cd(celldir) do
        jm = JointMap(;kwargs...)
        rf = get_response_fields(SpatialResponseFields, 10_000;kwargs...)
        jm, rf
    end
    res = conjunctions(jm, rf)
    FieldConjunctions{SpatialResponseFields}(rf, [res[k]["λ_infield"][:] for k in 1:length(res)],
                        [res[k]["λ_outfield"][:] for k in 1:length(res)],
                        [res[k]["sic_infield"] for k in 1:length(res)],
                        [res[k]["sic_outfield"] for k in 1:length(res)],
                        [res[k]["sic_sub"] for k in 1:length(res)])
end

function FieldConjunctions(::Type{T};do_save=true, redo=fname->false,kwargs...) where T <: AbstractResponseFields
end

function plot_conjunction(mm::SimpleMesh, place_field_idx, gaze_field_idx,res)
    m_floor = Translate(0.0, 0.0, -30.0)(floor_topology3())
    z_floor = zeros(nelements(m_floor))

    z_gaze = zeros(nelements(mm))
    z_gaze[gaze_field_idx] .= 1.0
    with_theme(plot_theme) do
        fig = Figure()
        ax1 = LScene(fig[2,1], show_axis=false)
        ax2 = LScene(fig[2,3], show_axis=false)
        ax3 = LScene(fig[2,5], show_axis=false)
        plotmesh!(ax1, mm;color=res["λ_outfield"], floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)
        fill!(z_floor, 1.0)
        z_floor[place_field_idx] .= 0.0
        viz!(ax1, m_floor;color=z_floor, colormap=:binary,showsegments=true, segmentcolor=:lightgray)
        Colorbar(fig[2,2], colorrange=extrema(filter(isfinite, res["λ_outfield"])), colormap=:binary,label="Firing rate [Hz]")
        plotmesh!(ax2, mm;color=res["λ_infield"], floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)
        fill!(z_floor, 0.0)
        z_floor[place_field_idx] .= 1.0
        viz!(ax2, m_floor;color=z_floor, colormap=:binary,showsegments=true, segmentcolor=:lightgray)
        Colorbar(fig[2,4], colorrange=extrema(filter(isfinite, res["λ_infield"])), colormap=:binary,label="Firing rate [Hz]")
        Label(fig[1,1], "Out of field", tellwidth=false)
        Label(fig[1,3], "In field", tellwidth=false)

        plotmesh!(ax3, mm;color=z_gaze, floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)


        link_cameras_lscene(fig)
        fig
    end
end

function plot_conjunction(fj::FieldConjunctions{SpatialResponseFields}, fields::GazeResponseFields, idx=1)
    nrefinements = fj.fields.args[:nrefinements]
    mm = get_maze_mesh(;nrefinements=nrefinements.g)
    m_floor = Translate(0.0, 0.0, -30.0)(floor_topology3(;nrefinements=nrefinements.p))
    z_floor = zeros(nelements(m_floor))
    clusters = merge_fields(fj.fields)
    place_field_idx = fj.fields.binidx[clusters[idx]]
    gaze_field_idx = fields.binidx
    z_gaze = zeros(nelements(mm))
    z_gaze[gaze_field_idx] .= 1.0
    with_theme(plot_theme) do
        fig = Figure()
        ax1 = LScene(fig[2,1], show_axis=false)
        ax2 = LScene(fig[2,3], show_axis=false)
        ax3 = LScene(fig[2,5], show_axis=false)
        plotmesh!(ax1, mm;color=fj.λ_outfield[idx], floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)
        fill!(z_floor, 1.0)
        z_floor[place_field_idx] .= 0.0
        viz!(ax1, m_floor;color=z_floor, colormap=:binary,showsegments=true, segmentcolor=:lightgray)
        Colorbar(fig[2,2], colorrange=extrema(filter(isfinite, fj.λ_outfield[idx])), colormap=:binary,label="Firing rate [Hz]")
        plotmesh!(ax2, mm;color=fj.λ_infield[idx], floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)
        fill!(z_floor, 0.0)
        z_floor[place_field_idx] .= 1.0
        viz!(ax2, m_floor;color=z_floor, colormap=:binary,showsegments=true, segmentcolor=:lightgray)
        Colorbar(fig[2,4], colorrange=extrema(filter(isfinite, fj.λ_infield[idx])), colormap=:binary,label="Firing rate [Hz]")
        Label(fig[1,1], "Out of field", tellwidth=false)
        Label(fig[1,3], "In field", tellwidth=false)
        Label(fig[1,5], "View and place fields", tellwidth=false)

        plotmesh!(ax3, mm;color=z_gaze, floor_offset=-15, ceiling_offset=10,colormap=:binary,showsegments=true,segmentcolor=:lightgray)
        viz!(ax3, m_floor;color=z_floor, colormap=:binary,showsegments=true, segmentcolor=:lightgray)

        link_cameras_lscene(fig)
        fig
    end
end