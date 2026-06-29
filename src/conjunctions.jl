# functions to compute conjunctions

struct FieldConjunctions{T1<:AbstractResponseFields, T2<:AbstractResponseFields}
    fields1::T1
    fields2::T2
    λ::Matrix{Float64}
    λ_infield::Matrix{Vector{Float64}}
    λ_outfield::Vector{Vector{Float64}}
end

function Hippocampus.issignificant(fj::FieldConjunctions;pv_threshold=0.05)
    res = fill(false, size(fj.λ_infield))
    for pidx in 1:length(fj.λ_outfield)
        for vidx in 1:size(fj.λ_infield,1)
            h = MannWhitneyUTest(fj.λ_outfield[pidx], fj.λ_infield[vidx,pidx])
            res[vidx,pidx] = pvalue(h;tail=:left) < pv_threshold
        end
    end
    res
end

function get_random_neighborhood(Am, i, k)
    aa = [i]
    while length(aa) < k
        idx = Am[:,i].nzind
        m = setdiff(idx,aa)
        if isempty(m)
            break
        end
        i = rand(m)
        push!(aa, i)
    end
    aa
end
#get_mesh(jf::FieldConjunctions{SpatialResponseFields}) = get_maze_mesh(;nrefinements=jf.fields.args[:nrefinements].g) 
#get_mesh(jf::FieldConjunctions{GazeResponseFields}) = Shadow("xy")(floor_topology3(;nrefinements=jf.fields.args[:nrefinements].p))

"""
Condition on either view or place fields and compare firing rates to the unconditioned

Compare the firing rates for e.g. view when the animal was in a particular place field vs 
    when it was elsewhere; if there is indeed a conjunction, the conditioned firing rate should be higher.
"""
function conjunctions(jm::JointMap, fields::SpatialResponseFields, viewidx::Vector{<:Integer})
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

function get_conditional_activity(jm::JointMap, fields1::T1, fields2::T2)  where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
    mm1 = get_mesh(T1, fields1.args[:nrefinements])
    mm2 = get_mesh(T2, fields2.args[:nrefinements])
    clusters1 = Hippocampus.merge_fields(fields1)
    nclusters1 = get_num_fields(fields1)
    clusters1 = clusters1[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    clusters2 = Hippocampus.merge_fields(fields2)
    nclusters2 = get_num_fields(fields2)
    clusters2 = clusters2[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]

    all_not_covered = setdiff(1:nelements(mm1), fields1.binidx)
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
    x_covered = zeros(nelements(mm2), length(clusters1))
    w_covered = zeros(nelements(mm2), length(clusters1))
    for (ii,covered) in enumerate(clusters1)
        # compare firing rates within vis outside the field
        n_covered = 0
        for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
            pidx = getindex(qidx,T1)
            vidx = getindex(qidx,T2)
            if pidx in covered
                #push!(Z_covered, λ)
                n_covered += 1
                x_covered[vidx,ii] += w
                w_covered[vidx,ii] += occ
            end
        end
    end
    (covered = (x=x_covered, w=w_covered), not_covered=(x=x_not_covered, w=w_not_covered))
end

"""
Return all valid indices om `idx` where the first index is not in `cluster1` but the second index is in `cluster2`
"""
function get_conditional_indices(idx::Vector{CartesianIndex{4}}, cluster1::Vector{<:Integer}, cluster2::Vector{<:Integer})
    valid_idx = Int64[]
    for (jj,ii) in enumerate(idx)
        k1 = getindex(ii,1)
        k2 = getindex(ii,2)
        if (k1 in cluster1) && (k2 in cluster2)
            push!(valid_idx,jj)
        end
    end
    valid_idx
end


function get_equivalent_cluster(idx::Vector{CartesianIndex{4}}, cluster1::Vector{<:Integer}, cluster2::Vector{<:Integer},n::Integer, D::AbstractMatrix{<:Real};dims=1)
    valid_idx = get_conditional_indices(idx, cluster1, cluster2)
    get_equivalent_cluster(idx,valid_idx, cluster1,cluster2,n,D;dims=dims)
end
"""
Create a cluster with `n` mmembers from from `cluster1` so that the second index is also in `cluster2`
"""
function get_equivalent_cluster(idx::Vector{CartesianIndex{4}}, valid_idx::Vector{<:Integer}, cluster1::Vector{<:Integer}, cluster2::Vector{<:Integer},n::Integer, D::AbstractMatrix{<:Real};dims=1,rng=Random.default_rng())
    clusters = [cluster1,cluster2]
    @assert length(clusters[dims]) == size(D,1)
    # grab a random index to use as a starting point
    fidx = rand(rng, valid_idx)
    # create a neighborhood
    if dims == 1
        midx = findfirst(cluster1.==getindex(idx[fidx],1))
        pidx = sortperm(D[:,midx])[1:n]
        new_cluster = cluster1[pidx]
    else
        midx = findfirst(cluster2.==getindex(idx[fidx],2))
        pidx = sortperm(D[:,midx])[1:n]
        new_cluster = cluster2[pidx]
    end
    new_cluster
end

function get_conditional_activity(jm::JointMap, cluster1::Vector{<:Integer}, cluster2::Vector{<:Integer},N::Union{Nothing, Int64}=nothing;dims=2)
    # get a set of valid indices, i.e. where elements of both cluster 1 and cluster 2 are present
    if dims==2
        dim_v = 1
        dim_p = 2
    else
        dim_v = 2
        dim_p = 1
    end
    if isnothing(N)
        N = maximum(getindex.(jm.index, dims))
    end
    X = zeros(N)
    W = zeros(N)
    f1 = in(cluster1)
    f2 = in(cluster2)
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        idx1 = getindex(qidx, 1)
        idx2 = getindex(qidx, 2)
        if f1(idx1)  && f2(idx2)
            if dims==2
                vidx = idx1
                pidx = idx2
            else
                vidx = idx2
                pidx = idx1
            end
            X[pidx] += w
            W[pidx] += occ
        end
    end
    X,W
end


# TODO: Create a more barebones version that is easier to test
function conjunctions_new(jm::JointMap, fields1::T1, fields2::T2;nshuffles=1000) where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
     mm1 = get_mesh(T1, fields1.args[:nrefinements])
    mm2 = get_mesh(T2, fields2.args[:nrefinements])
   
    D = distancematrix(mm1;between_centroids=false)
    # TODO: Get signifiance
    clusters1 = Hippocampus.merge_fields(fields1)
    nclusters1 = get_num_fields(fields1)
    clusters1 = clusters1[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    clusters2 = Hippocampus.merge_fields(fields2)
    nclusters2 = get_num_fields(fields2)
    clusters2 = clusters2[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]

    # distance matrix between points that are not in fields of fields1
    all_not_covered = setdiff(1:nelements(mm1), fields1.binidx)
    Da = D[all_not_covered, :]
    # compute outfield first
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

    λ_infield = zeros(nelements(mm2), length(clusters1))
    λ_covered = zeros(length(clusters2), length(clusters1))

    for (ii,covered_idx) in enumerate(clusters1)
        x_covered = zeros(nelements(mm2))
        w_covered = zeros(nelements(mm2))
        covered = fields1.binidx[covered_idx]
        for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
            pidx = getindex(qidx,T1)
            vidx = getindex(qidx,T2)
            if pidx in covered
                x_covered[vidx] += w
                w_covered[vidx] += occ
            end
        end
        λ_infield[:,ii] = x_covered./w_covered
        for (j,c) in enumerate(clusters2)
            qidx = fields2.binidx[c]
            _x = mean(x_covered[qidx])
            _w = mean(w_covered[qidx])
            λ_covered[j,ii] = _x./_w
        end
    end

    # now for the sub-sampling
    # grab indices for which element1 is the part of mm1 not covered by cluster 1
    # and element 2 in cluster 2
    
    λ_sub = fill(NaN, 1000, length(clusters2), length(clusters1))
    xx = fill(0.0, nelements(mm2))
    ww = fill(0.0, nelements(mm2))
    for (jj, _cidx2) in enumerate(clusters2)
        qidx = fields2.binidx[_cidx2]
        if T1 <: GazeResponseFields
            cidx1 = all_not_covered
            cidx2 = qidx 
            dims=1
        else
            cidx1 = qidx 
            cidx2 = all_not_covered
            dims=2
        end
        valid_idx = get_conditional_indices(jm.index, cidx1, cidx2)
        vidx1 = getindex.(jm.index[valid_idx],dims)
        # repeatedly grab an an index from 
        for (ii,_cidx1) in enumerate(clusters1)
            for kk in 1:nshuffles
                _vidx = rand(vidx1)
                # build a neighbourhood
                nnidx = all_not_covered[sortperm(Da[:,_vidx])[1:length(_cidx1)]]
                fill!(xx, 0.0)
                fill!(ww, 0.0)
                for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
                    pidx = getindex(qidx,T1)
                    vidx = getindex(qidx,T2)
                    if pidx in nnidx 
                        xx[vidx] += w
                        ww[vidx] += occ
                    end
                end
                λ_sub[kk,jj,ii] = mean(xx[qidx])/mean(ww[qidx])
            end
        end
    end
    λ_covered,λ_sub, λ_infield, λ_outfield
end

"""
    conjunctions2(jm::JointMap, viewidx::AbstractVector{<:Integer}, placeidx::AbstractVector{<:Integer})

Compute the mean firing rate when location was within `placeidx` and the view was in `viewidx`.

Compare to surrogates where the view was not in `view idx` by seelcting random views elicited by the same place idx, but from
outside the view field.
"""
# TODO: Test this!
function conjunctions2(jm::JointMap, viewidx::AbstractVector{<:Integer}, placeidx::AbstractVector{<:Integer}, condition_on::Integer=1)
    # get all view indices experienecd with placeidx
    matchindex(k1,k2) = condition_on == 1 ? (getindex(k1,2)==getindex(k2,2)) : (getindex(k1,1)==getindex(k2,1))
    ffp = in(placeidx)
    ffv = in(viewidx)
    ffc(idx) = ((condition_on==1 && ffv(getindex(idx,1))) || (condition_on==2 && ffp(getindex(idx,2))))
    ffk(idx) = ((condition_on==1 && ffp(getindex(idx,2))) || (condition_on==2 && ffv(getindex(idx,1))))
    in_view_idx = Int64[]
    out_view_idx = Int64[]
    time_in_field = 0.0
    spikes_in_field = 0.0
    time_out_field = 0.0
    spikes_out_field =0.0
    x_covered = 0.0
    w_covered = 0.0
    in_field_matched = Tuple{Int64, Int64, Float64, Float64}[]
    out_field_matched = Tuple{Int64, Int64, Float64, Float64}[]
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        vidx = getindex(qidx,1)
        pidx = getindex(qidx,2)
        if ffk(qidx)
        #if ffp(pidx)
            if ffc(qidx)
                push!(in_view_idx,vidx)
                time_in_field += occ
                spikes_in_field += w
                push!(in_field_matched, (vidx, pidx, w, occ))
            else
                # TODO: If there is a strong imbalance in gaze direction between
                # infield and outfield, we might get induced biases here
                push!(out_view_idx,vidx)
                time_out_field += occ
                spikes_out_field += w
                push!(out_field_matched, (vidx, pidx, w,occ))
            end
        end
    end
    λ_sub = fill(0.0, 1000)
    # TODO: Do we also need to worry about directional bias here?
    for r in 1:1000
        xx = 0.0
        ww = 0.0
        # draw one matched (v,p) from the outfield for each infield
        for vp in in_field_matched
            # find all outfield combos where the place idx was p
            midx = findall(k->matchindex(k,vp), out_field_matched)
            # grab a random index
            # what if we have no matches?
            if !isempty(midx)
                ridx = rand(midx)
                xx += out_field_matched[ridx][3]
                ww += out_field_matched[ridx][4]
            end
        end
        λ_sub[r] = xx/ww
    end
    λ_covered = spikes_in_field/time_in_field
    λ_covered, λ_sub
end

function get_spatial_rate_map(jm::JointMap, view_idx::AbstractVector{<:Integer},N::Integer)
    xx = zeros(N)
    yy = zeros(N)
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        vidx = getindex(qidx,1)
        pidx = getindex(qidx, 2)
        if vidx in view_idx
            xx[pidx] += w
            yy[pidx] += occ
        end
    end
    xx./yy
end

function get_view_rate_map(jm::JointMap, spatial_idx::AbstractVector{<:Integer},N::Integer)
    xx = zeros(N)
    yy = zeros(N)
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        vidx = getindex(qidx,1)
        pidx = getindex(qidx, 2)
        if pidx in spatial_idx
            xx[vidx] += w
            yy[vidx] += occ
        end
    end
    xx./yy
end

function conjunctions2(jm::JointMap, rf_gaze::GazeResponseFields, rf_spatial::SpatialResponseFields, condition_on::Integer)
    view_clusters = Hippocampus.merge_fields(rf_gaze)
    nclusters = Hippocampus.get_num_fields(rf_gaze)
    cidx = findall(dropdims(mean(nclusters,dims=2),dims=2).< 0.001)
    view_clusters = view_clusters[cidx]
    spatial_clusters = Hippocampus.merge_fields(rf_spatial)
    nclusters = Hippocampus.get_num_fields(rf_spatial)
    cidx = findall(dropdims(mean(nclusters,dims=2),dims=2).< 0.001)
    spatial_clusters = spatial_clusters[cidx]
    # compute infield/outfield
    if condition_on == 1
        # compute spatial fields conditioned on each of the view clusters
        N = length(rf_spatial.λ)
        all_outfield = setdiff(1:N, rf_spatial.binidx)
        λ_outfield = get_spatial_rate_map(jm, all_outfield,N)
        λ_infield = zeros(N, length(view_clusters))
        for (ii,vc) in enumerate(view_clusters)
            λ_infield[:,ii] = get_spatial_rate_map(jm, rf_gaze.binidx[vc],N)
        end
    else   
        N = length(rf_gaze.λ)
        all_outfield = setdiff(1:N, rf_gaze.binidx)
        λ_outfield = get_view_rate_map(jm, all_outfield,N)
        λ_infield = zeros(N, length(spatial_clusters))
        for (ii,sc) in enumerate(spatial_clusters)
            λ_infield[:,ii] = get_view_rate_map(jm, rf_spatial.binidx[sc],N)
        end
    end
    λ_covered = zeros(length(spatial_clusters), length(view_clusters))
    λ_sub = zeros(size(λ_covered)..., 1000)
    for (j,vc) in enumerate(view_clusters)
        for (i,sc) in enumerate(spatial_clusters)
            λ_covered[i,j],λ_sub[i,j,:] =Hippocampus.conjunctions2(jm, rf_gaze.binidx[vc], rf_spatial.binidx[sc],condition_on);
        end
    end
    λ_covered, λ_sub, λ_infield, λ_outfield
end

"""
Condition responses on `field1`
"""
function conjunctions(jm::JointMap, fields1::T1, fields2::T2;nshuffles=1000) where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
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
    for (ii,covered_idx) in enumerate(clusters1)
        covered = fields1.binidx[covered_idx]
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
        for (j,c) in enumerate(clusters2)
            qidx = fields2.binidx[c]
            _x = mean(x_covered[qidx])
            _w = mean(w_covered[qidx])
            λ_covered[j,ii] = _x./_w
        end
        for (jj,cluster) in enumerate(clusters2)
            qidx = fields2.binidx[cluster]
            for kk in 1:1000
                _x = mean(x_not_covered_sub[qidx,kk])
                _w = mean(w_not_covered_sub[qidx,kk])
                ll = _x/_w
                λ_sub[kk,jj,ii] = ll
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
    fname = "place_view_conjunction_new.jld2"
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
    nshuffles = get(kwargs, :nshuffles, 10_000)
    if do_compute
        jm = JointMap(;kwargs...)
        rf_spatial = get_response_fields(SpatialResponseFields,nshuffles;kwargs...)
        if isa(rf_spatial, JLD2.ReconstructedMutable)
            rf_spatial = get_response_fields(SpatialResponseFields,nshuffles;redo=fname->true, kwargs...)
        end
        rf_gaze = get_response_fields(GazeResponseFields,nshuffles;kwargs...)
        if isa(rf_gaze, JLD2.ReconstructedMutable)
            rf_gaze = get_response_fields(GazeResponseFields,nshuffles;redo=fname->true, kwargs...)
        end
        λ_covered, λ_sub,λ_infield, λ_outfield = conjunctions2(jm, rf_gaze, rf_spatial, 1)
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

function analyse(jm::JointMap, pvc::PlaceViewConjunction)
    FieldConjunctions(jm, pvc.spatial_fields, pvc.view_fields)
end

function marginalize(::Type{SpatialResponseFields}, X::Matrix{<:Real}, idx)
    dropdims(sum(X[:,idx],dims=2),dims=2)
end

function marginalize(::Type{GazeResponseFields}, X::Matrix{<:Real}, idx)
    dropdims(sum(X[idx,:],dims=1),dims=1)
end

function FieldConjunctions(jm::JointMap, fields1::T1, fields2::T2;smooth=true, smoothing_method=:laplace, α=0.1, niter=100,kwargs...) where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
    # TODO: Make sure we use the same bin filtering as for fields1 and fields2 here
    mm1 = get_mesh(T1, fields1.args[:nrefinements])
    mm2 = get_mesh(T2, fields2.args[:nrefinements])
    view_clusters = merge_fields(fields2)
    nclusters2 = get_num_fields(fields2)
    view_clusters = view_clusters[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]
    spatial_clusters = merge_fields(fields1)
    nclusters1 = get_num_fields(fields1)
    spatial_clusters = spatial_clusters[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    weight,occupancy = get_maps(jm)
    occupancy_v = dropdims(sum(occupancy,dims=2),dims=2)
    ps = dropdims(sum(occupancy,dims=1),dims=1)
    px = dropdims(sum(weight,dims=1),dims=1)
    # compute the (smoothed) firing rate in the original view field when conditoning on the place
    nvidx = setdiff(findall(occupancy_v.>0), fields2.binidx)
    λ_infield = Matrix{Vector{Float64}}(undef, length(view_clusters), length(spatial_clusters))
    λ_outfield = Vector{Vector{Float64}}(undef, length(spatial_clusters))
    λ = zeros(nelements(mm2), length(spatial_clusters))
    Ls = get_normalize_laplacian(mm2) 
    for (ii,sc) in enumerate(spatial_clusters)
        pidx = fields1.binidx[sc]
        # weighted average, since not all slices are weighted equally
        #xx = weight[:,pidx]*px[pidx]./sum(px[pidx])
        #xx = dropdims(sum(weight[:,pidx],dims=2),dims=2)
        xx = marginalize(T1, weight, pidx)
        ww = marginalize(T1, occupancy, pidx)
        #ww = dropdims(sum(occupancy[:,pidx],dims=2),dims=2)
        #ww = occupancy[:,pidx]*ps[pidx]./sum(ps[pidx])
        # smooth
        xs = laplace_smoothing(xx,Ls, α;niter=niter)
        ws = laplace_smoothing(ww,Ls, α;niter=niter)
        λ[:,ii] = xs./ws
        # Filter out unvisited bins< i.e. where the view occupancy is zero
        λ[occupancy_v.==0] .= NaN
        λ_outfield[ii] = xs[nvidx]./ws[nvidx]
        for (jj,vc) in enumerate(view_clusters)
            vidx = fields2.binidx[vc]
            λ_infield[jj,ii] = xs[vidx]./ws[vidx]
        end
    end
    FieldConjunctions{T1,T2}(fields1, fields2, λ, λ_infield, λ_outfield)
end

function FieldConjunctions(::Type{T1},::Type{T2};do_save=true, redo=fname->false,kwargs...) where T1 <: AbstractResponseFields where T2 <: AbstractResponseFields
    jm = JointMap(;redo=fname->false, do_save=true, kwargs...)
    rf1 = get_response_fields(T1 ,get(kwargs, :nshuffles, 10_000);kwargs...)
    rf2 = get_response_fields(T2 ,get(kwargs, :nshuffles, 10_000);kwargs...)
    FieldConjunctions(jm, rf1, rf2;kwargs...)
end

# Type for doing conditional information; we need this so that we can easily managed the shuffling
struct PlaceAccountingView{T<:Real}
    weight::Matrix{T}
    occupancy::Matrix{T}
    ee::T
    ees::Vector{T}
    args::Dict{Symbol,Any}
end

function process_kwargs(PlaceAccountingView, h::UInt32=zero(UInt32);nshuffles=1000, kwargs...)
    h = process_kwargs(JointMap,h;kwargs...)
    h = crc32c(string(:nshuffles=>nshuffles),h)
    h
end

function PlaceAccountingView(vpvrp::ViewAndPlaceRepresentationNew, jocc::JointOccupancy, jocc_filtered::JointFilteredOccupancy;nshuffles=1000, kwargs...)
    args = Dict{Symbol,Any}()
    for (k,v) in kwargs
        args[k] = v
    end
    args[:snuffles] = nshuffles
    jm = JointMap(vpvrp, jocc, jocc_filtered)
    method = get(kwargs, :smoothing_method, :laplace)
    α = get(kwargs, :α, 0.1)
    niter = get(kwargs, :niter, 50)
    jml = JointSmoothedMap(jm;method=method, α=α,niter=niter);
    # compute information about place given view
    ee = get_conditional_information(jml.weight, jml.occupancy)
    ees = zeros(nshuffles)
    prog = Progress(nshuffles)
    for i in 1:nshuffles
        jms = JointMap(vpvrp, jocc, jocc_filtered;shuffle_view=true);
        jmsl = JointSmoothedMap(jms;method=method, α=α, niter=niter);
        ees[i] = get_conditional_information(jmsl.weight, jmsl.occupancy)
        ProgressMeter.next!(prog)
    end
    PlaceAccountingView(jml.weight, jml.occupancy, ee, ees, args) 
end

function PlaceAccountingView(;redo=fname=>false, do_save=true, kwargs...)
    fname = "place_accounting_view.jld2"
    h = process_kwargs(PlaceAccountingView;kwargs...)
    if h > 0
        hs = string(h,base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo(fname) && isfile(fname)
        obj = load_jld2(PlaceAccountingView, fname)
    else
        sessiondir = DPHT.get_level_path("session")
        qdata, jocc = cd(sessiondir) do
            qdata = UnityRaytraceData(;kwargs...)
            jocc = JointOccupancy(;kwargs...)
            qdata, jocc
        end
        jocc_filtered = JointFilteredOccupancy(jocc,qdata;kwargs...)
        vpvrp = ViewAndPlaceRepresentationNew(;kwargs...)
        obj = get_conditional_information(PlaceAccountingView, vpvrp, jocc, jocc_filtered;kwargs...)
        if do_save
            save_jld2(obj, fname)
        end
    end
    obj
end

struct ViewAccountingSpace
end

## plots
function plot_conjunction(pvc::PlaceViewConjunction,idx=1;smooth=true,smoothing_method=:laplace, α=0.1, niter=100)

    # find the significant clusters
    view_clusters = merge_fields(pvc.view_fields)
    nclusters1 = get_num_fields(pvc.view_fields)
    view_clusters = view_clusters[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    spatial_clusters = merge_fields(pvc.spatial_fields)
    nclusters2 = get_num_fields(pvc.spatial_fields)
    spatial_clusters = spatial_clusters[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]
    mm = get_mesh(GazeResponseFields, pvc.view_fields.args[:nrefinements])
    mm2 = floor_topology3(;nrefinements=pvc.spatial_fields.args[:nrefinements].p)
    # translate down
    mm2 = Translate(0.0, 0.0, -20.0)(mm2)
    bbc = find_boundary(mm2, pvc.spatial_fields.binidx[spatial_clusters[idx]])
    if smooth
        # this is a bit clunky; we need to recompute weight and occupancy separately
        jm = cd(pvc.spatial_fields.args[:dir]) do
            # TODO: Make sure we use the correct parameters here
            JointMap()
        end
        jml = JointSmoothedMap(jm;method=:laplace, α=0.1, niter=100)
        ps = dropdims(sum(jml.occupancy,dims=1),dims=1)
        px = dropdims(sum(jml.weight,dims=1),dims=1)
        ccidx = pvc.spatial_fields.binidx[spatial_clusters[idx]] 
        X = jml.weight[:,ccidx]*px[ccidx]
        Y = jml.occupancy[:,ccidx]*ps[ccidx]
         λ_infield = X./Y
        _,occupancy = get_maps(jm)
        occ =  occupancy[:,ccidx]*ps[ccidx]
        λ_infield[occ.==0] .= NaN

        uccidx = setdiff(1:nelements(mm2), pvc.spatial_fields.binidx)
        X = jml.weight[:,uccidx]*px[uccidx]
        Y = jml.occupancy[:,uccidx]*ps[uccidx]
         λ_outfield = X./Y
        _,occupancy = get_maps(jm)
        occ =  occupancy[:,uccidx]*ps[uccidx]
        λ_outfield[occ.==0] .= NaN

    else
        λ_infield = pvc.λ_infield[:,idx]
        λ_outfield = pvc.λ_outfield
    end
     A = issignificant(pvc, pv_threshold=0.01)[:,idx]
    pidx = findall(A)
    scolor = fill(:black, length(A))
    scolor[pidx] .= :red
    cr = extrema([filter(isfinite, λ_infield);filter(isfinite, λ_outfield)])
    with_theme(plot_theme) do
        fig = Figure()
        Label(fig[1,1], "In field", tellwidth=false)
        Label(fig[1,2], "Out of field", tellwidth=false)
        lscene1 = LScene(fig[2,1], show_axis=false)
        plotmesh!(lscene1, mm;color=λ_infield,colorrange=cr,showsegments=true)
        lscene2 = LScene(fig[2,2], show_axis=false)
        plotmesh!(lscene2, mm;color=λ_outfield,colorrange=cr)
        # indicate the original view fields
        for lscene in [lscene1, lscene2]
            for (kk,vc) in enumerate(view_clusters)
                bb = find_boundary(mm, pvc.view_fields.binidx[vc])
                viz!(lscene, bb;color=scolor[kk])
            end
            viz!(lscene, mm2;color=:lightgray)
            viz!(lscene, bbc;color=:black)
        end
        # indicate the place field
        link_cameras_lscene(fig)
        # separate axis to show distribution of firing rate within each field
        ax3 = Axis(fig[2,3])
        yy = Float64[]
        xx = Float64[]
        for ii in 1:size(pvc.λ_sub,2)
            _yy = filter(isfinite, pvc.λ_sub[:,ii,idx])
            append!(yy, _yy)
            append!(xx, fill(ii, length(_yy)))
        end
        boxplot!(ax3, xx,yy;show_outliers=false)
        scatter!(ax3, 1:size(pvc.λ_covered,1), pvc.λ_covered[:,idx], color=scolor)
        ax3.yaxisposition = :right
        ax3.leftspinevisible = false
        ax3.rightspinevisible = true
        ax3.ylabel = "Firing rate [Hz]"
        colsize!(fig.layout, 3, Relative(0.25))
        fig
    end
    #plot_conjunction(mm, pvc.spatial_fields.binidx[spatial_clusters[1]], pvc.view_fields.binidx[view_clusters[1]])
end

"""
Plot view conditioned on place
"""
function plot_conjunction(pvc::PlaceViewConjunction,idx=1;smooth=true,smoothing_method=:laplace, α=0.1, niter=100,floor_offset=-20,colormap=:rain,_plot_theme=plot_theme, mazecolor=:darkgray,kwargs...)
    # find the significant clusters
    view_clusters = merge_fields(pvc.view_fields)
    nclusters1 = get_num_fields(pvc.view_fields)
    view_clusters = view_clusters[dropdims(mean(nclusters1,dims=2),dims=2) .< 0.001]
    spatial_clusters = merge_fields(pvc.spatial_fields)
    nclusters2 = get_num_fields(pvc.spatial_fields)
    spatial_clusters = spatial_clusters[dropdims(mean(nclusters2,dims=2),dims=2) .< 0.001]
    mm2 = get_mesh(GazeResponseFields, pvc.view_fields.args[:nrefinements])
    bbc = find_boundary(mm2, pvc.view_fields.binidx[view_clusters[idx]])
    ppc = pvc.view_fields.binidx[view_clusters[idx]]
    # get all points not part of a view field
    nppc = reduce(vcat, [pvc.view_fields.binidx[view_clusters[c]] for c in 1:length(view_clusters)])
    nppc = setdiff(1:nelements(mm2), nppc)

    mm = floor_topology3(;nrefinements=pvc.spatial_fields.args[:nrefinements].p)
    # translate down
    mm = Translate(0.0, 0.0, floor_offset)(mm)
    if smooth
        # this is a bit clunky; we need to recompute weight and occupancy separately
        jm = cd(pvc.spatial_fields.args[:dir]) do
            # TODO: Make sure we use the correct parameters here
            JointMap()
        end
        # joint smooth
        # this might be OK just for visualizing 
        jml = JointSmoothedMap(jm;method=:laplace, α=0.1, niter=100)
        ps = dropdims(sum(jml.occupancy,dims=2),dims=2)
        px = dropdims(sum(jml.weight,dims=2),dims=2)
        ccidx = pvc.view_fields.binidx[view_clusters[idx]] 
        #X = dropdims(mean(jml.weight[ccidx,:],dims=1),dims=1)
        X = jml.weight[ccidx,:]'*px[ccidx]
        Y = jml.occupancy[ccidx,:]'*ps[ccidx]
        #Y = dropdims(mean(jml.occupancy[ccidx,:],dims=1),dims=1)
        #res = get_conditional_activity(jm, pvc.view_fields, pvc.spatial_fields)
        λ_infield = X./Y
        _,occupancy = get_maps(jm)
        occ =  occupancy[ccidx,:]'*ps[ccidx]
        λ_infield[occ.==0] .= NaN
        #Xuc = vec(laplace_smoothing(res.not_covered.x, Ls, α;niter=niter))
        #Yuc = vec(laplace_smoothing(res.not_covered.w, Ls, α;niter=niter))
        uccidx = setdiff(1:nelements(mm2), pvc.view_fields.binidx)
        #X = dropdims(mean(jml.weight[uccidx,:],dims=1),dims=1)
        #Y = dropdims(mean(jml.occupancy[uccidx,:],dims=1),dims=1)
        X = jml.weight[uccidx,:]'*px[uccidx]
        Y = jml.occupancy[uccidx,:]'*ps[uccidx]
        occ =  occupancy[uccidx,:]'*ps[uccidx]
        λ_outfield = X./Y
        λ_outfield[occ.==0] .= NaN
    else
        λ_infield = pvc.λ_infield[:,idx]
        λ_outfield = pvc.λ_outfield
    end
    A = issignificant(pvc, pv_threshold=get(kwargs, :pv_threshold, 0.01))[:,idx]
    pidx = findall(A)
    @show pidx
    scolor = fill(:gray55, length(A))
    scolor[pidx] .= [:orange, :firebrick, :salmon, :goldenrod2][1:length(pidx)]
    cr = extrema([filter(isfinite, λ_infield);filter(isfinite, λ_outfield)])
    with_theme(_plot_theme) do
        fig = Figure()
        Label(fig[1,1], "In field", tellwidth=false)
        Label(fig[1,2], "Out of field", tellwidth=false)
        lscene1 = LScene(fig[2,1], show_axis=false)
        viz!(lscene1, mm;color=mazecolor)
        viz!(lscene1, mm;color=λ_infield,colorrange=cr,colormap=colormap)
        lscene2 = LScene(fig[2,2], show_axis=false)
        viz!(lscene2, mm;color=mazecolor)
        viz!(lscene2, mm;color=λ_outfield,colorrange=cr, colormap=colormap)
        # indicate the original view fields
        Z_infield = fill(NaN, nelements(mm2))
        Z_infield[ppc] .= 1.0
        Z_outfield = fill(NaN, nelements(mm2))
        Z_outfield[nppc] .= 1.0
        for lscene in [lscene1, lscene2]
            for (kk,bc) in enumerate(spatial_clusters)
                bb = find_boundary(mm, pvc.spatial_fields.binidx[bc])
                viz!(lscene, bb;color=scolor[kk])
            end
        end
        for (lscene, Zq,cc) in zip([lscene1, lscene2],[Z_infield, Z_outfield],[:goldenrod1,:steelblue4])
            # indicate the view vields
            plotmesh!(lscene, mm2;alpha=0, showsegments=true, segmentcolor=mazecolor,floor_offset=-10, ceiling_offset=10)
            #viz!(lscene, bbc;color=:black)
            plotmesh!(lscene, mm2;color=Zq,showsegments=false, floor_offset=-10, ceiling_offset=10, colormap=[cc])
        end
        # indicate the place field
        link_cameras_lscene(fig)
        # separate axis to show distribution of firing rate within each field
        ax3 = Axis(fig[2,3])
        yy = Float64[]
        xx = Float64[]
        for ii in 1:size(pvc.λ_sub,2)
            _yy = filter(isfinite, pvc.λ_sub[:,ii,idx])
            append!(yy, _yy)
            append!(xx, fill(ii, length(_yy)))
        end
        boxplot!(ax3, xx,yy;show_outliers=false,color=:steelblue4)
        ax3.xticklabelsvisible = false
        ax3.xticksvisible = false
        ax3.bottomspinevisible = false
        scatter!(ax3, 1:size(pvc.λ_covered,1), pvc.λ_covered[:,idx], color=scolor)
        ax3.yaxisposition = :right
        ax3.leftspinevisible = false
        ax3.rightspinevisible = true
        ax3.ylabel = "Firing rate [Hz]"
        colsize!(fig.layout, 3, 100)
        fig
    end
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