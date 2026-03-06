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

function get_rate_map(spf::SpatialResponseFields,jm::JointMap)
    nrefinements = spf.args[:nrefinements]
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    mm = SpatialMapNew(jm, m_floor)
    if spf.args[:smooth]
        # TODO: Make this more general
        sml = SmoothedMap(mm;method=spf.args[:smoothing_method],α=spf.args[:α], niter=spf.args[:niter])
        return get_rate_map(sml)
    end
    get_rate_map(mm)
end

function get_rate_map(spf::GazeResponseFields,jm::JointMap)
    nrefinements = spf.args[:nrefinements]
    mh = get_maze_mesh(;nrefinements=nrefinements.g)
    mm = ViewMapNew(jm, mh)
    if spf.args[:smooth]
        # TODO: Make this more general
        sml = SmoothedMap(mm;method=spf.args[:smoothing_method],α=spf.args[:α], niter=spf.args[:niter])
        return get_rate_map(sml)
    end
    get_rate_map(mm)
end

Base.getindex(qidx::CartesianIndex{4}, ::Type{SpatialResponseFields}) = getindex(qidx,2)
Base.getindex(qidx::CartesianIndex{4}, ::Type{GazeResponseFields}) = getindex(qidx,1)

function conjunctions(jm::JointMap, fields1::T1, fields2::T2) where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
    # T1 <: SpatialResponseFields, T2 <: GazeResponseFields
    mm1 = get_mesh(T1, fields1.args[:nrefinements])
    mm2 = get_mesh(T2, fields2.args[:nrefinements])

    # TODO: Get signifiance
    clusters1 = Hippocampus.merge_fields(fields1)
    nclusters1 = get_num_fields(fields1)
    clusters1 = clusters1[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    clusters2 = Hippocampus.merge_fields(fields2)
    nclusters2 = get_num_fields(fields2)
    clusters2 = clusters2[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]

    all_not_covered = setdiff(1:nelements(mm1), fields1.binidx)
    res = Dict()
    λ_covered = zeros(length(clusters2), length(clusters1))
    λ_sub = zeros(1000, length(clusters2), length(clusters1))
    λ_infield = zeros(nelements(mm2), length(clusters1))
    x_not_covered = zeros(nelements(mm2))
    w_not_covered = zeros(nelements(mm2))
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
       pidx = getindex(qidx,T1)
       vidx = getindex(qidx,T2) 
       if pidx in all_not_covered
            x_not_covered[vidx] += w
            w_not_covered[vidx] += occ
       end
    end
    λ_outfield = x_not_covered./w_not_covered
    for (ii,covered) in enumerate(clusters1)
        not_covered = setdiff(1:nelements(mm1), covered)
        not_covered_sub = [shuffle(all_not_covered)[1:length(covered)] for _ in 1:1000]
        # compare firing rates within vis outside the field
        x_covered = zeros(nelements(mm2))
        w_covered = zeros(nelements(mm2))
        x_not_covered = zeros(nelements(mm2))
        w_not_covered = zeros(nelements(mm2))

        x_not_covered_sub = zeros(nelements(mm2),length(not_covered_sub))
        w_not_covered_sub = zeros(nelements(mm2), length(not_covered_sub))
        n_covered = 0
        for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
            pidx = getindex(qidx,T1)
            vidx = getindex(qidx,T2)
            if pidx in covered
                #push!(Z_covered, λ)
                n_covered += 1
                x_covered[vidx] += w
                w_covered[vidx] += occ
            elseif pidx in all_not_covered
                #push!(Z_not_covered, λ)
                for (jj,nc) in enumerate(not_covered_sub)
                    if pidx in nc
                        x_not_covered_sub[vidx,jj] += w
                        w_not_covered_sub[vidx,jj] += occ
                    end

                end
            end
        end
        λ_infield[:,ii] = x_covered./w_covered
        # now we have firing rates for each of view bins conditioned on a particular place field
        # aggregate within each of the view fields
        λ_covered[:,ii]  = [mean(filter(isfinite, x_covered[fields2.binidx[c]]./w_covered[fields2.binidx[c]])) for c in clusters2]
        for (jj,cluster) in enumerate(clusters2)
            qidx = fields2.binidx[cluster]
            for kk in 1:1000
                ll = x_not_covered_sub[qidx,kk]./w_not_covered_sub[qidx,kk]
                λ_sub[kk,jj,ii] = mean(filter(isfinite, ll))
            end
        end
    end
    λ_covered,λ_sub, λ_infield, λ_outfield
end

function process_kwargs(::Type{<:AbstractFieldConjunctions},h::UInt32=zero(UInt32);kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = process_kwargs(SpatialResponseFields,h;kwargs...)
    h = process_kwargs(GazeResponseFields,h;kwargs...)
    h
end

function PlaceViewConjunction(;redo=fname->false, do_save=true, load_only=false, kwargs...)
    fname = "place_view_conjunction.jld2"
    h = process_kwargs(PlaceViewConjunction;kwargs...)
    if h > 0
        hs = string(h,base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = true 
    if !redo(fname) && isfile(fname)
        X = load_jld2(PlaceViewConjunction,fname)
       if isa(X, JLD2.ReconstructedMutable)
            do_compute = true
        else
            do_compute = false
        end
    end
    if do_compute && load_only
        # we just want to skip here
        return nothing
    end
    if do_compute
        jm = JointMap(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields,10_000;kwargs...)
        if isa(rf_spatial, JLD2.ReconstructedMutable)
            rf_spatial = get_response_fields(SpatialResponseFields,10_000;redo=fname->true, kwargs...)
        end
        rf_gaze = get_response_fields(GazeResponseFields,10_000;kwargs...)
        if isa(rf_gaze, JLD2.ReconstructedMutable)
            rf_gaze = get_response_fields(GazeResponseFields,10_000;redo=fname->true, kwargs...)
        end
        λ_covered, λ_sub,λ_infield, λ_outfield = conjunctions(jm, rf_spatial, rf_gaze)
        X = PlaceViewConjunction(rf_spatial, rf_gaze, λ_covered, λ_sub,λ_infield, λ_outfield)
        if do_save
            save_jld2(X,fname)
        end
    end
    X
end

function ViewPlaceConjunction(;redo=fname->false, do_save=true, load_only=false, kwargs...)
    fname = "view_place_conjunction.jld2"
    h = process_kwargs(ViewPlaceConjunction;kwargs...)
    if h > 0
        hs = string(h,base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    do_compute = true 
    if !redo(fname) && isfile(fname)
        X = load_jld2(ViewPlaceConjunction,fname)
       if isa(X, JLD2.ReconstructedMutable)
            do_compute = true
        else
            do_compute = false
        end
    end
    if do_compute && load_only
        # we just want to skip here
        return nothing
    end
    if do_compute
        jm = JointMap(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields,10_000;kwargs...)
        if isa(rf_spatial, JLD2.ReconstructedMutable)
            rf_spatial = get_response_fields(SpatialResponseFields,10_000;redo=fname->true, kwargs...)
        end
        rf_gaze = get_response_fields(GazeResponseFields,10_000;kwargs...)
        if isa(rf_gaze, JLD2.ReconstructedMutable)
            rf_gaze = get_response_fields(GazeResponseFields,10_000;redo=fname->true, kwargs...)
        end
        λ_covered, λ_sub,λ_infield, λ_outfield = conjunctions(jm, rf_gaze, rf_spatial)
        X = ViewPlaceConjunction(rf_gaze, rf_spatial, λ_covered, λ_sub,λ_infield, λ_outfield)
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

function plot_conjunction(pvc::PlaceViewConjunction,idx=1)

    # find the significant clusters
    view_clusters = merge_fields(pvc.view_fields)
    nclusters1 = get_num_fields(pvc.view_fields)
    view_clusters = view_clusters[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    spatial_clusters = merge_fields(pvc.spatial_fields)
    nclusters2 = get_num_fields(pvc.spatial_fields)
    spatial_clusters = spatial_clusters[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]
    mm = get_mesh(GazeResponseFields, pvc.view_fields.args[:nrefinements])

    with_theme(plot_theme) do
        fig = Figure()
        Label(fig[1,1], "In field", tellwidth=false)
        Label(fig[1,2], "Out of field", tellwidth=false)
        lscene1 = LScene(fig[2,1], show_axis=false)
        plotmesh!(lscene1, mm;color=pvc.λ_infield[:,idx])
        lscene2 = LScene(fig[2,2], show_axis=false)
        plotmesh!(lscene2, mm;color=pvc.λ_outfield)
        # indicate the original view fields
        for lscene in [lscene1, lscene2]
            for vc in view_clusters
                bb = find_boundary(mm, pvc.view_fields.binidx[vc])
                viz!(lscene, bb;color=:black)
            end
        end
        fig
    end
    #plot_conjunction(mm, pvc.spatial_fields.binidx[spatial_clusters[1]], pvc.view_fields.binidx[view_clusters[1]])
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

function plot_conjunctions(pvc::PlaceViewConjunction)
    vm = get_rate_map(pvc.view_fields)
    spm = get_rate_map(pvc.spatial_fields)
end

function plot_conjunctions(fj::FieldConjunctions)
    (nrows,ncols) = size(fj.λ_infield)
    with_theme(plot_theme) do
        fig = Figure(size=(1267, 742))
        for c in 1:ncols
            for r in 1:nrows
                if nrows > ncols
                    lg = GridLayout(fig[c,r])
                else
                    lg = GridLayout(fig[r,c])
                end
                plot_conjunctions!(lg, fj, c,r)
            end
        end
        fig
    end

end

function plot_conjunctions(fj::FieldConjunctions,args...)
    with_theme(plot_theme) do
        fig = Figure()
        lg = GridLayout(fig[1,1])
        plot_conjunctions!(lg, fj,args...)
        fig
    end
end

function plot_conjunctions!(lg, fj::FieldConjunctions{SpatialResponseFields,GazeResponseFields},pidx=1,vidx=1;indicate_view_field=true)
    lscene = LScene(lg[1,1], show_axis=false)
    mm = get_mesh(GazeResponseFields,fj.fields2.args[:nrefinements])
    plotmesh!(lscene, mm;color=:lightblue, ceiling_offset=10, floor_offset=-10, colormap=:binary)
    plotmesh!(lscene, mm;color=fj.λ[:,pidx], ceiling_offset=10, floor_offset=-10, colormap=:binary)
    if indicate_view_field
        plot_response_fields!(lscene, fj.fields2, vidx;ceiling_offset=10, floor_offset=-10)
    end
    # indicate the place field we are using
    m_floor = Translate(0.0, 0.0, -25)(Hippocampus.floor_topology3(;nrefinements=fj.fields1.args[:nrefinements].p))
    clusters = merge_fields(fj.fields1)
    Z = fill(NaN, nelements(m_floor))
    Z[fj.fields1.binidx[clusters[pidx]]] .= 1.0
    viz!(m_floor;color=:lightblue)
    viz!(m_floor;color=Z,colormap=:binary)
    # show distribution of rates within and outside of the view fields
    ax = Axis(lg[1,2])
    ax.yaxisposition = :right
    ax.bottomspinevisible = false
    ax.xticksvisible = false
    ax.yticksvisible = true
    ax.xticks = ([1,2], ["Outfield", "Infield"])
    ax.xticklabelrotation = -π/5
    ax.leftspinevisible = false
    ax.rightspinevisible = true
    ax.ylabel = "Firing rate [Hz]"
    xx = fill(1.0, length(fj.λ_outfield[pidx]))
    append!(xx, fill(2.0, length(fj.λ_infield[vidx,pidx])))
    yy = fj.λ_outfield[pidx]
    append!(yy, fj.λ_infield[vidx,pidx])
    bbx = boxplot!(ax, xx, yy, color=:darkgray,show_outliers=false, show_notch=true)
    colsize!(lg, 2, 100)
    # TODO: Indicate significance
    h = MannWhitneyUTest(fj.λ_outfield[pidx], fj.λ_infield[vidx,pidx])
    pv = pvalue(h;tail=:left)
    if pv < 0.001
        tt = "**"
    elseif pv < 0.01
        tt = "*"
    else
        tt = "ns"
    end
    ymax = maximum(bbx.q5s[])
    ymin = minimum(bbx.q1s[])
    Δy = ymax-ymin
    ylims!(ax, ymin-0.1*Δy, ymax+0.2*Δy)
    bracket!(ax, 1, ymax, 2, ymax, offset=5, text=tt,style=:square)
end

function plot_conjunctions!(lg, fj::FieldConjunctions{GazeResponseFields,SpatialResponseFields},pidx=1,vidx=1;indicate_view_field=true)
    lscene = LScene(lg[1,1], show_axis=false)
    mm = get_mesh(GazeResponseFields,fj.fields1.args[:nrefinements])
    plotmesh!(lscene, mm;color=:lightblue, ceiling_offset=10, floor_offset=-10, colormap=:binary)
    
    # indicate the place field we are using
    m_floor = Translate(0.0, 0.0, -25)(Hippocampus.floor_topology3(;nrefinements=fj.fields2.args[:nrefinements].p))
    clusters = merge_fields(fj.fields1)
    Z = fill(NaN, nelements(mm))
    Z[fj.fields1.binidx[clusters[pidx]]] .= 1.0
    plotmesh!(lscene, mm;color=Z, ceiling_offset=10, floor_offset=-10, colormap=:binary)

    viz!(m_floor;color=:lightblue)
    viz!(m_floor;color=fj.λ[:,pidx],colormap=:binary)
    if indicate_view_field
        plot_response_fields!(lscene, fj.fields2, vidx;offset=-25)
    end
    # show distribution of rates within and outside of the view fields
    ax = Axis(lg[1,2])
    ax.yaxisposition = :right
    ax.bottomspinevisible = false
    ax.xticksvisible = false
    ax.yticksvisible = true
    ax.xticks = ([1,2], ["Outfield", "Infield"])
    ax.xticklabelrotation = -π/5
    ax.leftspinevisible = false
    ax.rightspinevisible = true
    ax.ylabel = "Firing rate [Hz]"
    xx = fill(1.0, length(fj.λ_outfield[pidx]))
    append!(xx, fill(2.0, length(fj.λ_infield[vidx,pidx])))
    yy = fj.λ_outfield[pidx]
    append!(yy, fj.λ_infield[vidx,pidx])
    bbx = boxplot!(ax, xx, yy, color=:darkgray,show_outliers=false, show_notch=true)
    colsize!(lg, 2, 100)
    # TODO: Indicate significance
    h = MannWhitneyUTest(fj.λ_outfield[pidx], fj.λ_infield[vidx,pidx])
    pv = pvalue(h;tail=:left)
    if pv < 0.001
        tt = "**"
    elseif pv < 0.01
        tt = "*"
    else
        tt = "ns"
    end
    ymax = maximum(bbx.q5s[])
    ymin = minimum(bbx.q1s[])
    Δy = ymax-ymin
    ylims!(ax, ymin-0.1*Δy, ymax+0.2*Δy)
    bracket!(ax, 1, ymax, 2, ymax, offset=5, text=tt,style=:square)
end