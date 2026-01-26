# functions to compute conjunctions

"""
Condition on either view or place fields and compare firing rates to the unconditioned

Compare the firing rates for e.g. view when the animal was in a particular place field vs 
    when it was elsewhere; if there is indeed a conjunction, the conditioned firing rate should be higher.
"""
function conjunctions(jm::JointMap, fields, mm::SimpleMesh, m_floor::SimpleMesh)
    # TODO: What should the fields be? 
    # We actually need to condition
    covered = Set{Int64}()
    for field in fields
        bidx = findall(Meshes.intersects.(m_floor, field))
        for b in bidx
            push!(covered, b)
        end
    end
    not_covered = setdiff(1:nelements(mm), covered)
    # compare firing rates within vis outside the field
    x_covered = zeros(nelements(mm))
    w_covered = zeros(nelements(mm))
    x_not_covered = zeros(nelements(mm))
    w_not_covered = zeros(nelements(mm))
    for (w,occ,qidx) in zip(jm.weight, jm.occupancy, jm.index)
        pidx = getindex(qidx,2)
        vidx = getindex(qidx,1)
        if pidx in covered
            #push!(Z_covered, λ)
            x_covered[vidx] += w
            w_covered[vidx] += occ
        elseif pidx in not_covered
            #push!(Z_not_covered, λ)
            x_not_covered[vidx] += w
            w_not_covered[vidx] += occ
        end
    end
    λ_covered  = x_covered./w_covered
    λ_not_covered = x_not_covered./w_not_covered
    λ_covered, λ_not_covered
end