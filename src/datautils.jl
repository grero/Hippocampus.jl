mutable struct RunCounter
    level::Int64
end

function (lc::RunCounter)(fname::String)
    lm = lc.level > 0
    lc.level -= 1
    return lm > 0
end

struct LevelCounter
    origin::String
    idx::Int64
    depth::Int64
end

function LevelCounter(origin::String, d)
    idx = findfirst(DPHT.levels.==origin)
    LevelCounter(origin, idx, d)
end

function (lc::LevelCounter)(fname::String,cwd=pwd())
    ll = DPHT.level(cwd)
    idx = findfirst(DPHT.levels.==ll)
    d = lc.idx-idx
    return d <= lc.depth 
end

function process_multi_session_cell(celldir::String)
    if occursin("session", celldir)
        return celldir # this is not a multi-session
    end
    daydir = DPHT.get_level_path("day", celldir)
    unity_eyelinkfile = joinpath(daydir,"session01","unityfile_eyelink_new.csv")
    target_link = joinpath(daydir,"unityfile_eyelink_new.csv")
    if isfile(unity_eyelinkfile) && !islink(target_link)
        symlink(unity_eyelinkfile, target_link)
    end
    # This is a bit weird
    dayname = DPHT.get_level_name("session", celldir)
    celldirs = glob(replace(celldir, dayname=>joinpath(dayname, "session??")))
    combined_dir = joinpath(daydir, "session_combined")
    if !isdir(combined_dir)
        mkdir(combined_dir)
    end
    for _celldir in celldirs
        sn = DPHT.get_level_name("session", _celldir)
        targetlink = replace(_celldir, sn=>"session_combined")
        if !isdir(targetlink)
            mkpath(targetlink)
        end
    end
    new_celldir = replace(celldir, dayname=>joinpath(dayname, "session_combined"))
    @assert isdir(new_celldir)
    new_celldir
end