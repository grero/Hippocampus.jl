find_place_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_place_selective, celldirs;kwargs...)
find_view_selective_cells(celldirs::Vector{String};kwargs...) = find_selective_cells(is_view_selective, celldirs;kwargs...)

function find_selective_cells(func::Function, celldirs::Vector{String};skip_error=true, kwargs...)
    _is_selective = fill(false, length(celldirs))
    for (i,celldir) in enumerate(celldirs)
        try
            _is_selective[i] = func(celldir;kwargs...)
        catch ee
            if !skip_error
                rethrow(ee)
            end
        finally
            continue
        end
    end
    _is_selective
end

function is_place_selective(celldir::String;filtering="All", raycast_type="1px")
    fname=joinpath(celldir, "Filt$(filtering)",raycast_type,"vmpc.mat")
    if !ispath(fname)
        error("No vmpc object found at $fname")
    end

    vp = MAT.matread(fname)
    vp_data = vp["vmp"]["data"]
    is_selective(vp_data)
end

function is_view_selective(celldir::String;filtering="All", raycast_type="1px")
    fname=joinpath(celldir, "Filt$(filtering)",raycast_type,"vmsv.mat")
    if !ispath(fname)
        error("No vmsv object found at $fname")
    end

    vp = MAT.matread(fname)
    vp_data = vp["vms"]["data"]
    is_selective(vp_data)
end
    
function is_selective(vp_data::Dict{String,Any})
    # check for 
    if "SIC" in keys(vp_data)
        sic_value = vp_data["SIC"]
    elseif "SIC_adsm" in keys(vp_data)
        sic_value = vp_data["SIC_adsm"]
    else
        error("No SIC variable found")
    end
    if "SICsh" in keys(vp_data)
        sic_threshold = percentile(vec(vp_data["SICsh"]), 95)
    elseif "crit_sm" in keys(vp_data)
        sic_threshold = vp_data["crit_sm"]
    else
        error("No shuffled data found")
    end

    return sic_value > sic_threshold 
end