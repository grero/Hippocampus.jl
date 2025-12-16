using MAT
using CSV
using DataFrames
using Makie
using ProgressMeter

abstract type AbstractViewMap <: AbstractMap end
struct DummyCam
    pos::Point3f
    dir::Vec3f
    fov::Float32
    z_near::Float32
    frustrum_ratio::Float32
end

"""
Trace a ray from position `x,y` through the camera with focal length `focal_length` until it
impacts something in the arena
"""
function raytrace(x, y, cam::DummyCam,mm::MazeModel)
    # TODO: It looks like the raytracing function in Unity just uses the viewport. In other words,
    # what the camera 'sees' is a normalized coordinate system (not the physical sensor.)
    # find the angle of the point
    # height of frustrum at near clip
    # Unity is also using the far instead of the near plane, setting the far plane at 25 units
    fovr = cam.fov
    fwidth = tan(fovr/2)*cam.z_near
    fheight = fwidth/cam.frustrum_ratio
    # is the normalized plane from -1 to 1?
    xc = x*fwidth
    yc = y*fheight
    # flip x since x-values left-of-center should be associated with a positive angle
    θ = atan(xc, cam.z_near)
    θc = atan(cam.dir[2],cam.dir[1])
    θ += θc
    # TODO: Assumes no camera elevation angle
    ϕ = atan(yc,cam.z_near) 
    v = [cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ)]
    # alternatively
    #b = [cam.dir nullspace(permutedims(cam.dir))]
    #v = xc*b[:,2] + yc*b[:,3]
    #v = v + cam.dir
    # re-normalize
    #v = v./norm(v)
    
    dl = 0.001
    (xp,yp,zp) = cam.pos
    while true
        dx,dy,dz = dl*v 
        xp += dx
        yp += dy 
        zp += dz 
        if impacts([xp,yp,zp], mm)
            xp -= dx
            yp -= dy
            zp -= dz
            break
        end
    end
    xp,yp,zp
end

function projecto(cam::DummyCam, pos::AbstractVector{T}) where T <: Real
    fwidth = tan(cam.fov/2)*cam.z_near
    fheight = fwidth/cam.frustrum_ratio
    # project along camera axis
    # TODO: Not clear which is which of these axes
    b = [cam.dir nullspace(permutedims(cam.dir))]
    pos_p = b'*(pos-cam.pos)
    xh = (pos_p[2]*cam.z_near/pos_p[1])/fwidth
    yh = (pos_p[3]*cam.z_near/(sqrt(pos_p[1]^2+pos_p[2]^2)))/fheight
    xh,yh
end

struct GazeOnMaze
    time::Vector{Vector{Float64}}
    gaze::Vector{Matrix{Float64}}
    fixation::Vector{Vector{Bool}}
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    header::Dict
end

numtrials(gdata::GazeOnMaze) = length(gdata.gaze)

DPHT.filename(::Type{GazeOnMaze}) = "maze_raytrace.mat"
DPHT.level(::Type{GazeOnMaze}) = "session"

function get_trial(gdata::GazeOnMaze, i;trial_start=1)
    tg = gdata.time[i]
    gg = gdata.gaze[i]
    fm = gdata.fixation[i]
    tg,gg,fm
end

# TODO: Make sure that this actually works
function GazeOnMaze(;do_save=true, redo=false)
    fname = DPHT.filename(GazeOnMaze)
    if !redo && isfile(fname)
        qdata = MAT.matread(fname)
        args = Any[]
        for k in fieldnames(GazeOnMaze)
            push!(args, qdata[string(k)])
        end
        gdata = GazeOnMaze(args...)
    else
        edata = cd(DPHT.process_level(EyelinkData)) do
            EyelinkData()
        end
        udata = cd(DPHT.process_level(UnityData)) do
            UnityData()
        end
        gdata = GazeOnMaze(edata, udata)
        if do_save
            DPHT.save(gdata)
        end
    end
    gdata
end

function DPHT.save(gdata::T;append_tag=true) where T
    fname = DPHT.filename(T)
    qdata = Dict{String,Any}()
    metadata = Dict{String,Any}() 
    if append_tag
        tag!(metadata, storepatch=true)
    end
    for k in fieldnames(T)
        v = getfield(gdata, k)
        qdata[string(k)] = v
    end
    qdata["meta"] = metadata
    MAT.matwrite(fname, qdata)
end

function save_jld2(gdata::T,fname=DPHT.filename(T), ;append_tag=true) where T
    if SimpleMesh in fieldtypes(T)
        # this is a bit hacky, needed because SimpleMesh does not save cleanly
        qdata = Dict()
        for k in fieldnames(T)
            if fieldtype(T, k) <: SimpleMesh
                continue
            end
            v = getfield(gdata, k)
            qdata[k] = v
        end
    else
        qdata = gdata
    end
    fname = replace(fname, ".mat"=>".jld2")
    metadata = Dict{String,Any}() 
    if append_tag
        tag!(metadata, storepatch=true)
    end
    JLD2.save(fname, Dict("data"=>qdata, "meta"=>metadata))
end

function load_jld2(::Type{T},fname=DPHT.filename(T)) where T
    fname = replace(fname, ".mat"=>".jld2")
    meta,data = JLD2.load(fname, "meta","data")
    ft = fieldtypes(T)
    if SimpleMesh in ft 
        midx = findfirst(ft.==SimpleMesh)
        args = Any[]
        for k in fieldnames(T)
            if fieldtype(T, k) <: SimpleMesh
                continue
            end
            v = data[k]
            push!(args, v)
        end
        #hackis
        if T <: SpatialMapNew
            mm = Shadow("xy")(floor_topology3())
        else
            mm = get_maze_mesh()
        end
        insert!(args, midx, mm)
        mdata = T(args...)
    else
        mdata = data
    end
    mdata
end

function DPHT.load(::Type{T}) where T
    fname = DPHT.filename(T)
    fname_jld2 = replace(fname ,".mat"=>".jld2")
    if isfile(fname_jld2)
        return load_jld2(T)
    end
    qdata = MAT.matread(fname)
    metadata = qdata["meta"]
    args = Any[]
    for k in fieldnames(T)
        push!(args, qdata[string(k)])
    end
    T(args...)
end

"""
Convert `x` from pixel coordinates to sensor coordinates
"""
function scale_to_camera(x, sensor_width, screen_width)
    x = sensor_width*(x -0.5*screen_width)/screen_width
end

"""
Trace a ray from position `x,y` through the camera with focal length `focal_length` until it
impacts something in the arena
"""
function raytrace(x, y, pos,direction, fov, near_clip=0.3;camera_height=1.85,frustrum_ratio=1.78)
    # TODO: It looks like the raytracing function in Unity just uses the viewport. In other words,
    # what the camera 'sees' is a normalized coordinate system (not the physical sensor.)
    # find the angle of the point
    # height of frustrum at near clip
    fovr = π*fov/180
    fheight = 2*tan(fovr/2)*near_clip
    fwidth = fheight*frustrum_ratio
    # is the normalized plane from -1 to 1?
    xc = x*fwidth
    yc = y*fheight
    # flip x since x-values left-of-center should be associated with a positive angle
    θ = atan(xc, near_clip)
    ϕ = atan(yc,near_clip) 
    θ += direction
    # now cast along the line (θ,ϕ)

    dl = 0.01
    xp,yp,zp = (pos[1], pos[2],camera_height)
    sθ = sin(θ)
    cθ = cos(θ)
    sϕ = sin(ϕ)
    while true
        dx,dy,dz = (dl*sθ, dl*cθ,dl*sϕ)
        xp += dx
        yp += dy 
        zp += dz 
        if impacts([xp,yp,zp])
            xp -= dx
            yp -= dy
            zp -= dz
            break
        end
    end
    xp,yp,zp
end

"""
Construct GazeMaze from eyelinkata and raytrace data from Unity
"""
function GazeOnMaze(edata::EyelinkData, raytracedata::DataFrame)
    # get the trial triggers in raytrace time
    triggers = edata.triggers
    time_raytrace = raytracedata.Column2[:]
    gaze_pos = [raytracedata.Column10 raytracedata.Column12 raytracedata.Column11]
    gaze_pos[ismissing.(gaze_pos)] .= NaN
    gaze_pos = something.(gaze_pos)
    Δt = Float64(time_raytrace[1]) - edata.timestamps[1,1]
    new_timestamps = edata.timestamps .+ Δt
    nt = size(edata.triggers,1)
    gaze = Vector{Matrix{Float64}}(undef, nt)
    gtime = Vector{Vector{Float64}}(undef,nt) 
    fixation = Vector{Vector{Bool}}(undef, nt)
    for i in 1:nt
        idx0 = searchsortedfirst(time_raytrace, new_timestamps[i,1])
        idx1 = searchsortedlast(time_raytrace, new_timestamps[i,3])
        gaze[i] = permutedims(gaze_pos[idx0:idx1,:])
        # convert back to eyelink time
        gtime[i] = time_raytrace[idx0:idx1] .- time_raytrace[idx0] .+ edata.timestamps[1,1]
        gtime[i] = (gtime[i] .- gtime[i][1])/1000.0 # convert to seconds
        # TODO: Actually set this
        fixation[i] = fill(false, idx1-idx0+1)
    end
    GazeOnMaze(gtime,gaze,fixation,triggers,edata.timestamps,Dict())
end

"""
Return the 3D eye position on objects in the maze
"""
function GazeOnMaze(edata::EyelinkData, udata::UnityData)
    # we an only align edata and udata using events, so we need to operate on trials
    # the most compact representation
    nt = size(udata.triggers,1)
    gaze = Vector{Matrix{Float64}}(undef, nt)
    gtime = Vector{Vector{Float64}}(undef,nt)
    fixation = Vector{Vector{Bool}}(undef, nt)
    screen_width = edata.header["gaze_coords"][3] - edata.header["gaze_coords"][1]
    screen_height = edata.header["gaze_coords"][4] - edata.header["gaze_coords"][2]
    gx0, gy0, gxm, gym = edata.header["gaze_coords"]
    # find center of screen
    gx0 = (gxm-gx0)/2
    gy0 = (gym-gy0)/2
    prog = Progress(nt,desc="Raytracing...")
    for i in 1:nt
        # get the time index of the eyelink data triggers
        _t_e,gx,gy,fixation[i] = get_trial(edata, i)

        # Unity camera sensor size: (36,24)
        t_u,posx,posy,dir = get_trial(udata,i)

        # use the same reference
        # convert eyelink time to seconds
        t_e = (_t_e .- _t_e[1])/1000.0
        t_u .= t_u .- t_u[1]

        _gaze = fill(0.0, 3,length(t_e))
        gtime[i] = t_e
        for (j,t) in enumerate(t_e)
            tidx = searchsortedlast(t_u,t)
            outofbounds = false
            if gx[j] < 0.0 || gx[j] > screen_width
               _gaze[1,j] = NaN
               outofbounds = true
            end
            if gy[j] < 0.0 || gy[j] > screen_width
               _gaze[2,j] = NaN
               outofbounds=true
            end
            if !outofbounds
                #x = scale_to_camera(gx[j]-gx0, 36.0, screen_width)
                #y = scale_to_camera(gy[j]-gy0, 24.0, screen_height)
                # scale to viewport coordinates
                x = (gx[j] - gx0)/screen_width
                y = (gy[j] - gy0)/screen_height
                _gaze[:,j] .= raytrace(x,y,[posx[tidx],posy[tidx]],π*dir[tidx]/180,60.0, 0.3)
            end
        end
        gaze[i] = _gaze
        next!(prog)
    end
    GazeOnMaze(gtime,gaze,fixation, edata.triggers, edata.timestamps, Dict("focal_length"=>50.0,"camera_height"=>2.5))
end

function Makie.convert_arguments(::Type{<:AbstractPlot}, gdata::GazeOnMaze)
    # 3D scatter plot of all positions
    nn = sum([size(g,2) for g in gdata.gaze])
    gaze = Vector{Point3f}(undef, nn)
    offset = 0
    for (ii,gaze) in enumerate(gdata.gaze)
        n = size(gaze,2)
        gaze[offset+1:offset+n] = Point3f.(eachcol(gaze))
        offset += n
    end
    ax3 = S.Axis3(plots=[S.Scatter(gaze)])
    # TODO: Use the "exploded view here", that is indicate the pillars, as well as the floor
    # and ceiling
    S.GridLayout(ax3)
end

function visualize!(lscene, gdata::GazeOnMaze;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0),fixation_only=true, kwargs...)
    nt = numtrials(gdata)
    current_gaze = Observable([Point3f(NaN)])
    gdata_trial = lift(trial) do _trial
        if 0 < _trial.i <= nt
            tg = gdata.time[_trial.i]
            return tg, gdata.gaze[_trial.i],gdata.fixation[_trial.i]
        else
            return Float64[], fill(0.0, 3, 0), Bool[]
        end
    end

    current_j = 1
    onany(gdata_trial, current_time) do _gdt, _ct
        tg = _gdt[1] 
        gaze = _gdt[2]
        if fixation_only
            fixmask = _gdt[3]
        else
            fixmask = fill(true, length(tg))
        end
        j = searchsortedfirst(tg, _ct) 
        if 0 < j <= length(fixmask)
            _current_j = current_j
            current_j = j
            j0 = min(_current_j, j)
            j1 = max(_current_j, j)
            _fixmask = fixmask[j0:j1]
            current_gaze[] = Point3f.(eachcol(gaze[:,j0:j1][:,_fixmask]))
        else
            current_gaze[] = [Point3f(NaN)]
        end
    end

    scatter!(lscene, current_gaze, color=RGB(0.8, 0.8, 0.8))
end

struct UnityRaytraceData
    analogtime::Vector{UInt64}
    positions::Matrix{Float32}
    directions::Vector{Float32}
    triggers::Matrix{Union{UInt64,Missing}}

    gaze::Vector{Matrix{Float64}}
    position::Vector{Matrix{Float64}}
    head_direction::Vector{Vector{Float64}}
    timestamps::Vector{Vector{Float64}}
    fixated_object::Vector{Vector{String}}
    fixating::Vector{Vector{Bool}}
    unityfile::String
end

function UnityRaytraceData(analogtime, positions, direction, gaze::Vector{Vector{Matrix{Float64}}}, position, head_direction, timestamps,fixated_object,fixating)
    nt = length(gaze)
    triggers = fill(typemax(UInt64), nt, 3)
    UnityRaytraceData(analogtime, positions, direction, triggers, gaze, position, head_direction, timestamps, fixated_object, fixating, "unityfile_eyelink.csv")
end
"""
    get_gaze(X::UnityRaytraceData)

Return a matrix of all gaze positons concatenated
"""
function get_gaze(X::UnityRaytraceData;only_fixations=true)
    nn = div.(length.(X.gaze),3)
    nt = sum(nn)
    Y = zeros(Float64, 3, nt)
    w = zeros(Float64, nt)
    offset = 0
    for (_fix,_gaze,tt) in zip(X.fixating, X.gaze, X.timestamps)
        if only_fixations
            fidx = _fix
        else
            fidx = fill(true, size(_gaze,2))
        end
        Δt = diff(tt)
        push!(Δt, mean(Δt))
        nf = sum(fidx)
        Y[:,offset+1:offset+nf] .= _gaze[:,fidx]
        w[offset+1:offset+nf] .= Δt[fidx]
        offset += nf
    end
    Y[:,1:offset], w[1:offset]
end

function UnityRaytraceData(;do_save=true, redo=false,append_tag=true, raytrace_fname="unityfile_eyelink.csv", extradir::String="",fix_eyelink=false, kwargs...)
    fname = DPHT.filename(UnityRaytraceData)
    if !redo && isfile(fname)
        t1 = time()
        ut = DPHT.load(UnityRaytraceData)
        t2 = time()-t1
        @debug "Time" t2
        if !isfile(replace(fname, ".mat"=>".jld2"))
            save_jld2(ut)
        end
    else
        edata = cd(DPHT.process_level(EyelinkData)) do
            EyelinkData()
        end
        if !isempty(extradir)
            raytrace_fname = joinpath(extradir, raytrace_fname)
        end
        if !ispath(raytrace_fname)
            error("$(raytrace_fname) not found in $(pwd())")
        end
        unity_eyelinkfile = CSV.File(raytrace_fname, header=0)
        n = length(unity_eyelinkfile)
        fixated_points = fill(NaN, 3, n)
        position = fill(NaN32, 3, n)
        direction = fill(NaN32, n)
        timestamps = zeros(UInt64,n)
        fixated_object = fill("unknown", n)
        i = 1
        for row in unity_eyelinkfile
            # TODO: Grab more data here
            timestamps[i] = row[2]
            px,py,pz,a = (row[6],row[7],row[8],row[9])
            if any(ismissing.((px,py,pz,a)))
                continue
            end
            θ = π*a/180.0
            # FIXME: This does not appear to be the actual eyelink timestamp!
            gx,gy,gz = (row[10],row[11],row[12])
            if (gx !== missing) && (gy !== missing) && (gz !== missing)
                # gx,gy,gz is relative to player ?
                # swap z and y
                position[:,i] .= (px,pz,py)
                direction[i] = θ
                fixated_object[i] = row[3]

                fixated_points[:,i] .= (gx,gz,gy) # unity has the z-axis into the scene
                i += 1
            end
        end

        timestamps .-= timestamps[1]

        # break up into trials using edata
        nt = numtrials(edata)
        trial_fixations = Vector{Matrix{Float64}}(undef,nt)
        trial_position = Vector{Matrix{Float64}}(undef, nt)
        trial_head_direction = Vector{Vector{Float64}}(undef, nt)
        trial_times = Vector{Vector{Float64}}(undef,nt)
        trial_fixated_object = Vector{Vector{String}}(undef, nt)
        fixating = Vector{Vector{Bool}}(undef, nt)
        #t0 = edata.analogtime[1]
        te,_,_, = get_trial(edata, 1)
        t0 = te[1]
        triggers = edata.timestamps .- t0
        for i in 1:nt
            if any(ismissing.(edata.triggers[i,:]))
                trial_fixations[i] = zeros(3,0)
                trial_position[i] = zeros(3,0) 
                trial_head_direction[i] = zeros(0)
                trial_times[i] = zeros(0)
                trial_fixated_object[i] = String[]
                fixating[i] = Bool[]
                continue
            end
            te,_,_,fm = get_trial(edata, i;trial_start=2)
            # unit raytraced data use time relative to start of recording
            te .-= t0
            idx0 = searchsortedfirst(timestamps, te[1])
            idx1 = searchsortedlast(timestamps,te[end])
            trial_fixations[i] = fixated_points[:,idx0:idx1]
            trial_position[i] = position[:,idx0:idx1]
            trial_head_direction[i] = direction[idx0:idx1]
            trial_times[i] = timestamps[idx0:idx1]/1000.0 # convert to seconds
            trial_fixated_object[i] = fixated_object[idx0:idx1]
            # we need to match fixations to the actual time points
            _fm = fill(false, idx1-idx0+1)
            for j in 1:length(_fm)
                _idx = searchsortedlast(te, timestamps[idx0+j-1])
                if 0 < _idx < length(te) 
                    _fm[j] = fm[_idx]
                end
            end
            fixating[i] = _fm
        end
        ut = UnityRaytraceData(timestamps, position, direction,triggers, trial_fixations, trial_position, trial_head_direction, trial_times, trial_fixated_object, fixating, raytrace_fname)
        if do_save
            save_jld2(ut;append_tag=append_tag)
            #DPHT.save(ut;append_tag=append_tag)
        end
    end
    return ut
end

numtrials(gdata::UnityRaytraceData) = length(gdata.gaze)

DPHT.filename(::Type{UnityRaytraceData}) = "unity_raytrace.mat"
DPHT.level(::Type{UnityRaytraceData}) = "session"

function visualize!(lscene, unitygaze::UnityRaytraceData;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0),indicate_object=false,kwargs...)
    ugdata_trial = lift(trial) do _trial
        #Point3f.(eachcol(unitygaze.gaze[_trial.i]))
        gaze = unitygaze.gaze[_trial.i]
        pos = unitygaze.position[_trial.i]
        dir = unitygaze.head_direction[_trial.i]
        tu = unitygaze.timestamps[_trial.i]
        tu .-= tu[1]
        tg = [Point3f(gaze[[1,2,3],i]) for i in 1:size(gaze,2)]
        tp = [Point3f(pos[[1,2,3],i]) for i in 1:size(gaze,2)]
        td = dir
        tu,tg,tp,td, unitygaze.fixated_object[_trial.i] 
    end

    current_pos = Observable(ugdata_trial[][3][1:1])
    current_arrow_pos = Observable(ugdata_trial[][3][1:1])
    current_gaze = Observable(ugdata_trial[][2])
    current_path = Observable(ugdata_trial[][3])
    current_dir = Observable(ugdata_trial[][4][1:1])
    fixated_object = Observable(ugdata_trial[][5][1])
    current_arrow = lift(current_dir) do θ
        [mean([Point3f(sin(_θ), cos(_θ), 0.0) for _θ in θ])]
    end
    current_j = 1

    onany(ugdata_trial, current_time) do _ugt, _ct
        _tu, gpoints,ppoints,_td = (_ugt[1], _ugt[2], _ugt[3],_ugt[4])
        j = searchsortedfirst(_tu, _ct) 
        if j <= length(_tu)
            j0 = min(j,current_j)
            j1 = max(j, current_j)
            current_pos[] = ppoints[j0:j1]
            current_dir[] = _td[j0:j1]
            current_arrow_pos[] = [mean(ppoints[j0:j1])]
            current_gaze[] = gpoints[j0:j1]
            current_path[] = ppoints
            fixated_object[] = mode(_ugt[5][j0:j1])
            current_j = j
        end
    end
    lines!(lscene, current_path)
    scatter!(lscene,current_pos, color=:black)
    arrows!(lscene, current_arrow_pos, current_arrow, color=:black)
    scatter!(lscene, current_gaze, color=:red)
    if indicate_object
        text!(lscene, 0.5, 0.85;text=fixated_object,space=:relative)
    end
end

"""
Wrapper type to visualize the raytracing
"""
abstract type AbstractRaytraceViewer end;

struct RaytraceViewer <: AbstractRaytraceViewer
    udata::UnityData
    gazemaze::GazeOnMaze
end

function get_trial(raytrace::RaytraceViewer)
    tu,px,py,hh = get_trial(raytrace.udata, _trial.i)
    tg, ga = get_trial(raytrace.gazemaze, _trial.i)
    tg,ga,tu,px,py,hh
end

struct UnityRaytraceViewer <: AbstractRaytraceViewer
    raytrace_data::UnityRaytraceData
end

function get_trial(raytrace::UnityRaytraceViewer,i::Integer;trial_start=1)
    tt = raytrace.timestamps[i]
    tg = raytrace.gaze[i]
    fixmask = fill(true, size(tg,2))
    tt,tg,fixmask
end

function get_trial(raytrace::UnityRaytraceData,i::Integer;trial_start=1,not_on=["HintImage","CueImage","Robot"])
    tt = raytrace.timestamps[i]
    tg = raytrace.gaze[i]
    tp = raytrace.position[i]
    hd = raytrace.head_direction[i]
    fixmask = raytrace.fixating[i]
    fo = raytrace.fixated_object[i]
    # filter out objects that we don't want to include, amending fixmask accordingly
    fixmask .&= (!in(not_on)).(raytrace.fixated_object[i])
    tt,tg,tp,fixmask,fo,hd
end

function visualize!(lscene, raytrace::RaytraceViewer;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0),kwargs...)
    data_trial = lift(trial) do _trial
        tu,px,py,hh = get_trial(raytrace.udata, _trial.i)
        tg, ga = get_trial(raytrace.gazemaze, _trial.i)
        tg,ga,tu,px,py,hh
    end
    current_ray = Observable([Point3f(NaN)=>Point3f(NaN)])
    ray_color = Observable(parse(Colorant, :green))
    current_jg = 1
    onany(data_trial, current_time) do _ugt, _ct
        tg,ga,tu,px,py,hh = _ugt
        tu .-= tu[1]
        tg .-= tg[1]
        ju = searchsortedfirst(tu, _ct) 
        if 0 < ju <= length(tu)
            jg = searchsortedfirst(tg, _ct)
            j0 = max(min(current_jg, jg),0)
            j1 = min(max(current_jg, jg),length(tg))
            current_jg = jg
            current_ray[] = [Point3f(px[ju],py[ju],1.85)=>Point3f(ga[:,j]) for j in j0:j1]
            # whether the gaze vector is within the field of view
            xx = ga[1,j0:j1] .- px[ju]
            yy = ga[2,j0:j1] .- py[ju]
            zz = ga[3,j0:j1] .- 1.85
            θ = atan.(yy,xx)
            ϕ = atan.(zz,xx./sin.(θ))
            if any((ϕ .< -π/6 .|| ϕ .> ϕ./6) .&& (cos.(θ .- hh[ju]) .> cos(π/6)))
                ray_color[] = parse(Colorant, :red)
            else
                ray_color[] = parse(Colorant, :green)
            end
        end
    end
    linesegments!(lscene, current_ray, color=ray_color)
end

function compute_histogram(gdata::GazeOnMaze,mm::MazeModel;fixations_only=true)
    bins = get_bins(mm)
    if fixations_only 
        gaze = Vector{Matrix{Float64}}(undef, length(gdata.gaze))
        weight = Vector{Vector{Float64}}(undef, length(gdata.gaze))
        for i in eachindex(gaze)
            Δt = diff(gdata.time[i])
            fix = gdata.fixation[i]
            fix[end] = false
            weight[i] = Δt[fix[1:end-1]]
            gaze[i] = gdata.gaze[i][:,fix]
        end
    else
        gaze = Vector{Matrix{Float64}}(undef, length(gdata.gaze))
        weight = Vector{Vector{Floaft64}}(undef, length(gdata.gaze))
        for i in eachindex(gaze)
            Δt = diff(gdata.time[i])
            weight[i] = Δt
            gaze[i] = gdata.gaze[i][:,1:end-1]
        end
    end
    counts = Dict{Symbol,Vector{Array{Float64,3}}}()
    for k in keys(bins)
        counts[k] = compute_histogram(gaze,bins[k],weight)
    end
    counts,bins
end

function compute_histogram(gdata::UnityRaytraceData,mm::MazeModel)
    bins = get_bins(mm)
    gaze = Vector{Matrix{Float64}}(undef, length(gdata.gaze))
    weight = Vector{Vector{Float64}}(undef, length(gdata.gaze))
    for i in eachindex(gaze)
        Δt = diff(gdata.timestamps[i])
        weight[i] = Δt
        gaze[i] = gdata.gaze[i][:,1:end-1] # Skip the last point
    end

    counts = Dict{Symbol,Vector{Array{Float64,3}}}()
    idx = Vector{Vector{Tuple{Int64,Int64,Int64,Symbol}}}(undef, length(gdata.gaze))
    # initialize
    for ii in 1:length(idx)
        idx[ii] = fill((0,0,0,:unknown), size(gaze[ii],2))
    end
    for k in keys(bins)
        counts[k] = compute_histogram!(idx, gaze,bins[k],weight,k)
    end
    # replace the object place holder with the actual fixated object as reported
    # by unity
    for ii in 1:length(idx)
        fo = gdata.fixated_object[ii]
        _idx = idx[ii]
        for jj in 1:length(idx[ii])
            _idx[jj] = (_idx[jj][1:3]..., Symbol(fo[jj]))
        end
    end
    counts,bins,idx
end

"""
Type to indicate that we want to replay the maze as the subject experienced it.
"""
struct MazeReplayer
    mm::MazeModel
    udata::UnityData
    camera_aspect::Float32
end

MazeReplayer(mm::MazeModel, udata::UnityData) = MazeReplayer(mm,udata,1.0f0)

function visualize!(lscene::LScene, mp::MazeReplayer;current_time::Observable{Float64}=Observable(0.0), trial::Observable{Trial}=Observable(Trial(1)), kwargs...)

    visualize!(lscene,mp.mm;show_ceiling=true)
    aspect = mp.camera_aspect 
    cam = Makie.camera(lscene.scene)
    cam.projection[] = Makie.perspectiveprojection(60.0f0, aspect, 0.3f0, 1000.0f0)
    cc = cameracontrols(lscene.scene)
    #cc.settings.clipping_mode[] = :adaptive # :static
    udata = mp.udata
    nt = numtrials(udata)
    udata_trial = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return tp,posx,posy,dir = get_trial(udata,_trial.i)
        else
            return Float64[], Float64[], Float64[], Float[]
        end
    end

    onany(udata_trial, current_time) do _udt, _ct
        tp = _udt[1]
        tp .-= tp[1]
        j = searchsortedfirst(tp, _ct)
        px,py,dir = _udt[2:end]
        if 0 < j <= length(tp)
            pos = Point3f(px[j],py[j], 1.85)
            θ = π*dir[j]/180.0 # convert to radians
            cc.lookat[] = Point3f(0.3*sin(θ), 0.3*cos(θ), 0.0) + pos
            cc.eyeposition[] = pos
            update_cam!(lscene.scene, cc)
        end
    end
end

#TODO: Add "follow" mode
function show_maze(args...;show_axis=false, kwargs...)
    fig = Figure()
    lscene = LScene(fig[1,1], show_axis=show_axis)
    show_maze!(lscene, args...;kwargs...)
    fig
end

function show_maze!(lscene, bins,counts::Union{Dict{Symbol,Vector{Array{T,3}}},Nothing}=nothing,normals::Union{Dict{Symbol,Vector{Vector{Float64}}},Nothing}=nothing;explore=false, replay=false, interactive=false, gdata::Union{Nothing, GazeOnMaze}=nothing, udata::Union{Nothing, UnityData}=nothing, trial::Observable{Int64}=Observable(1),trialtime::Observable{Float64} = Observable(0.0), offsets::Union{Nothing, Dict{Symbol, Vector{Vector{Float64}}}}=nothing,show_ceiling=true,posters=nothing) where T <: Real
    #ax = Axis3(fig[1,1],aspect=:data)
    for k in keys(bins)
        for (i,bin) in enumerate(bins[k])
            if k == :ceiling && show_ceiling == false
                continue
            end
            m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
            if offsets !== nothing && k in keys(offsets)
                m = Translate(offsets[k][i]...)(m)
            end
            if counts !== nothing
                n = normals[k][i]
                # we want to color only the inside
                c = counts[k][i]
                _color = fill!(similar(c), 0.0)
                for d in 1:length(n)
                    if n[d] < 0
                        idx = ntuple(dim->dim==d ? 1 : axes(c,dim), 3)
                        _color[idx...] .= dropdims(sum(c,dims=d),dims=d)
                    elseif n[d] > 0
                        idx = ntuple(dim->dim==d ? size(c,d) : axes(c,dim), 3)
                        _color[idx...] .= dropdims(sum(c,dims=d),dims=d)
                    end
                end
                _color = _color[:]
            else
                if k in keys(pillar_color)
                    _color = pillar_color[k]
                else
                    _color = RGB(0.8, 0.8, 0.8) 
                end
            end
            viz!(lscene, m, color=_color,colormap=:Blues)
        end
    end
    if posters !== nothing
        # use the positions from the udata header
        if udata !== nothing
            _poster_pos = udata.poster_pos
        else
            _poster_pos = poster_pos 
        end
        wall_idx = assign_posters(bins, normals)
        rot = LinearMap(RotX(3π/2))   
        for (ii,(pp,img)) in enumerate(zip(_poster_pos,posters))
            sp = sprite(img, Rect2(-1.25, -2.5/1.2/2, 2.5, 2.5/1.2))
            sp2 = rot(sp)
            trans = LinearMap(Translation(pp[1],pp[2], 2.5))
            θ = acos(sp2.normals[1]'*normals[wall_idx[ii][1]][wall_idx[ii][2]])
            rot2 = LinearMap(RotZ(θ))
            sp3 = trans(rot2(sp2))
            plot!(lscene, sp3)
        end
    end

    udata_trial = lift(trial) do _trial
        if udata !== nothing
            tp,posx,posy,dir = get_trial(udata,_trial)
            return tp .- tp[1], permutedims([posx posy]),dir
        else
            return Float64[], fill(0.0, 3, 0), Float64[]
        end
    end

    gdata_trial = lift(trial) do _trial
        if gdata !== nothing
            return gdata.time[_trial], gdata.gaze[_trial],gdata.fixation[_trial]
        else
            return Float64[], fill(0.0, 3, 0), Bool[]
        end
    end
    
    position = lift(udata_trial) do _udata
        [Point3f(pos[1], pos[2], 0.5) for pos in eachcol(_udata[2])]
    end

    ii = Observable(1)

    gaze_pos = Observable([Point3f(NaN)])
    current_j = 1

    current_gaze = lift(ii) do i
        tp = udata_trial[][1]
        tg = gdata_trial[][1] 
        gaze = gdata_trial[][2]
        fixmask = gdata_trial[][3]
        j = searchsortedfirst(tg, tp[i]) 
        _fixmask = fixmask[current_j:j]
        _current_j = current_j
        current_j = j
        Point3f.(eachcol(gaze[:,_current_j:j][:,_fixmask]))
    end

    
    lookat = Point3f(1.0, 0.0, 2.5)
    on(events(lscene.scene).scroll, priority=20) do (dx,dy)
        i_new = round(Int64,ii[] + 5*dx)
        if 0 < i_new <= length(position[])
            ii[] = i_new
        end
    end

    if replay
        # replay experiment with the supplied position
        #cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        cc = cameracontrols(lscene.scene)
        cc.fov[] = 60.0 
        if gdata !== nothing
            tg,gaze,fixmask = (gdata.time[trial[]], gdata.gaze[trial[]],gdata.fixation[trial[]])
            tg .-= tg[1]
        end

        #tp,posx,posy,head_direction = get_trial(udata,trial[])
        #tp .-= tp[1]
        #position = permutedims([posx posy])
        

        on(ii) do i
            # grab the points 
            tp, position, head_direction = udata_trial[]
            _tp = tp[i]
            pos = Point3f(position[1,i], position[2,i], 2.5)
            θ = π*head_direction[i]/180

            # only use fixation points
            #gaze_pos[] = [Point3f(dropdims(mean(gaze[:,current_j+1:j][:,_fixmask],dims=2),dims=2))]
            if gdata !== nothing
                _fixmask = fixmask[current_j:j]
                j = searchsortedfirst(tg, _tp)
                gaze_pos[] = Point3f.(eachcol(gaze[:,current_j:j][:,_fixmask]))
                current_j = j
            end
            cc.lookat[] = Point3f(cos(θ), sin(θ), 0.0) + pos
            cc.eyeposition[] = pos
            update_cam!(lscene.scene, cc)
        end
        scatter!(lscene, gaze_pos, color=:red)

        if !interactive
            @async for j in 1:length(udata_trial[][1])
                ii[] = j
                yield()
                sleep(0.03)
            end
        else
            ii[] = 1 
        end
    elseif explore
        # set up camera inside of the maze
        cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        #cc.eyeposition[] = Point3f(0.0, 0.0, 2.5)
        eyepos = Point3f(0.0, 0.0, 2.5)
        v = lookat - eyepos 
        v = v./norm(v)
        #translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0,2.5))
        update_cam!(lscene.scene, eyepos, lookat)
        on(events(lscene.scene).keyboardbutton, priority=20) do event
            if ispressed(lscene.scene, Keyboard.up)
                pos = cc.eyeposition[]
                dx = 0.1*v
                npos = pos + dx
                translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0, -0.1))
                if impacts(cc.eyeposition[])
                    # move back
                    # TODO: This doesn't quite work, but maybe we don't care
                    translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0, 0.1))
                    # last coordinate if foward movement (for some inexplicable reason))
                end
            end
            if ispressed(lscene.scene, Keyboard.right)
                rotate_cam!(lscene.scene, cc, Point3f(0.0, -0.1, 0.0))
                return Consume()
            end
            if ispressed(lscene.scene, Keyboard.left)
                rotate_cam!(lscene.scene,cc, Point3f(0.0, 0.1, 0.0))
                return Consume()
            end
        end
    else
        #if gdata !== nothing
        #   tg,gaze,fixmask = (gdata.time[trial[]], gdata.gaze[trial[]],gdata.fixation[trial[]])
        #   scatter!(lscene, gaze[1,fixmask], gaze[2,fixmask], gaze[3,fixmask],color=:red)
        #end
        
        # show current position
        current_pos = lift(ii) do i
            if udata !== nothing
                pos = position[][i]
            else
                pos = Point3f(NaN)
            end
            [pos]
        end

        current_dir = lift(ii) do i
            head_direction = udata_trial[][3]
            if udata !== nothing
                θ = π*head_direction[i]/180
                return [Point3f(cos(θ), sin(θ), 0.0)]
            end
            return [Poin3f(NaN)]
        end
       
        lines!(lscene, position,color=:black)
        scatter!(lscene, current_pos, color=:blue)
        scatter!(lscene, current_gaze, color=:black)
        arrows!(lscene, current_pos, current_dir,color=:blue)
    end
end

function explore_maze(mm::MazeModelNew, points::Vector{Point3f}=Point3f[])
    with_theme(plot_theme) do
        fig = Figure()
        lscene = LScene(fig[1,1])
        plot!(lscene, mm)
        # add lights
        set_lights!(lscene, [])
        push_light!(lscene, PointLight(RGBf(1,1,1), Point3f(-9.1, 5.0, 1.5)))
        push_light!(lscene, PointLight(RGBf(1,1,1), Point3f(-4.8, -9.0, 1.5)))
        push_light!(lscene, PointLight(RGBf(1,1,1), Point3f(5.0, 9.0, 1.5)))
        push_light!(lscene, PointLight(RGBf(1,1,1), Point3f(9.0, -5.0, 1.5)))
        push_light!(lscene, PointLight(RGBf(1,1,1), Point3f(5.0, 1.1, 1.5)))
        if !isempty(points)
            scatter!(lscene, points)
        end

        #set up camera
        lookat = Point3f(1.0, 0.0, 0.7)
        cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        eyepos = Point3f(0.0, 0.0, 0.7)
        v = lookat - eyepos 
        v = v./norm(v)
        #translate_cam!(lscene.scene, cc, Point3f(0.0, 0.0,2.5))
        update_cam!(lscene.scene, eyepos, lookat)
        on(events(lscene.scene).keyboardbutton, priority=20) do event

            if ispressed(lscene.scene, Keyboard.up) || ispressed(lscene.scene, Keyboard.down)
                if ispressed(lscene.scene, Keyboard.up)
                    dx = Point3f(0.0, 0.0, -0.1)
                else
                    dx = Point3f(0.0, 0.0, 0.1)
                end
                translate_cam!(lscene.scene, cc, dx)
                if impacts(cc.eyeposition[], mm)
                    # move back
                    # TODO: This doesn't quite work, but maybe we don't care
                    translate_cam!(lscene.scene, cc, -dx)
                    # last coordinate if foward movement (for some inexplicable reason))
                end
            end
            if ispressed(lscene.scene, Keyboard.right)
                rotate_cam!(lscene.scene, cc, Point3f(0.0, -0.1, 0.0))
                return Consume()
            end
            if ispressed(lscene.scene, Keyboard.left)
                rotate_cam!(lscene.scene,cc, Point3f(0.0, 0.1, 0.0))
                return Consume()
            end
        end
        fig
    end
end

# these are the types we can visualize
Visualizables = Union{MazeModel, UnityData}

# TODO: It would be more elegant to make use of Makie recipe here

# I feel like this is duplicating functionality that must be in Makie somwhere...
function create_axis(::Type{T},fig;kwargs...) where T <: Makie.AbstractAxis
    if T <: LScene
        axis_args = (show_axis=get(kwargs, :show_axis,false),)
    else
        axis_args = (backgroundcolor=get(kwargs, :backgroundcolor, :white),)
    end
    lscene = T(fig[1,1];axis_args...)
end


function create_axis(obj::EyelinkData, fig;kwargs...)
    axtype = get_axis_type(EyelinkData)
    ax = create_axis(axtype, fig;kwargs...)
    # hide everything
    hidedecorations!(ax)
    ax.backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.0) 
    ax
end

function create_axis(obj::UnityData, fig;kwargs...)
    axtype = get_axis_type(UnityData)
    ax = create_axis(axtype, fig;kwargs...)
    ax
end

function create_axis(obj::MazeReplayer, fig;kwargs...)
    axtype = get_axis_type(MazeReplayer)
    ax = create_axis(axtype, fig;kwargs...)
end

function create_axis(obj::MazeModel, fig;kwargs...)
    axtype = get_axis_type(MazeModel)
    ax = create_axis(axtype, fig;kwargs...)
end

function create_axis(obj::Posters, fig;kwargs...)
    axtype = get_axis_type(Posters)
    ax = create_axis(axtype,fig;kwargs...)
end

 get_axis_type(::Type{T}) where T <: Any = LScene

function visualize(objects;kwargs...)
    fig = Figure()
    
    # attach events
    current_time = Observable(0.0)
    on(events(fig.scene).scroll, priority=20) do (dx,dy)
        current_time[] = current_time[] + dx
    end
    current_trial = Observable(Trial(1))
    on(events(fig.scene).keyboardbutton, priority=20) do event
        has_changed = false
        nc = current_trial[].i
        if ispressed(fig.scene, Keyboard.up)
            nc += 1
            # TODO: Is there a meaninful way to check whether we have reached the end here?
            has_changed = true
        elseif ispressed(fig.scene, Keyboard.down)
            nc  = current_trial[].i
            if nc-1 > 0
                has_changed = true
                nc -= 1
            end
        end
        if has_changed
            current_time[] = 0.0
            current_trial[] = Trial(nc)
        end
    end
    lg = GridLayout(fig[1,1])
    # show title with current trial and current time
    stitle = lift(current_trial, current_time) do ctrial,ctime
        "Trial: $(ctrial.i) time: $ctime"
    end
    Label(lg[1,1],stitle, tellheight=true, tellwidth=false)
    visualize(lg[2,1], objects;current_time=current_time, trial=current_trial, kwargs...)
    current_trial[] = Trial(1)
    current_time[] = 0.0
    fig
end

VectorOrMatrix{T} = Union{Vector{T}, Matrix{T}}

function visualize(fig::Union{Figure,GridLayout,GridPosition}, objects::VectorOrMatrix{T};kvs...) where T <: Union{T2, NTuple{N,T2},Vector{T2}} where T2 where N
    # tuple indicates overlay, i.e. plot in the same axis
    # vector means stack, i.e. axis occupying the same grid locaiton
    scenes = Any[]
    lg = GridLayout(fig[1,1])
    fobjects = Any[]
    scene_offset = 0
    for jj in axes(objects,2)
        for ii in axes(objects,1)
            obj = objects[ii,jj]
            if isa(obj,Vector)
                for (kk,_obj) in obj
                    axtype = get_axis_type(typeof(_obj))
                    push!(scenes, create_axis(axtype, lg[ii,jj];kvs...))
                    push!(fobjects, _obj)
                end
                scene_offset += length(obj)
            else
                if isa(obj, Tuple)
                    axtypes = Any[]
                    _obj = first(obj)
                    axtype = get_axis_type(typeof(_obj))
                    _scene = create_axis(_obj, lg[ii,jj];kvs...)
                    push!(axtypes, axtype)
                    push!(fobjects, _obj)
                    push!(scenes,_scene)
                    for _obj in obj[2:end]
                        axtype = get_axis_type(typeof(_obj))
                        kk = findfirst(axtypes.==axtype)
                        if kk !== nothing
                            #FIXME: This is wrong, scenes contains all the scenes, not just for this cell
                            _scene = scenes[scene_offset+kk]
                            push!(scenes, _scene)
                        else
                            _scene = create_axis(_obj, lg[ii,jj];kvs...)
                            push!(axtypes, axtype)
                            push!(scenes, _scene)
                        end
                        push!(fobjects, _obj)
                    end
                    scene_offset += length(obj)
                else
                    axtype = get_axis_type(typeof(obj))
                    push!(scenes, create_axis(obj, lg[ii,jj];kvs...))
                    scene_offset += 1
                    push!(fobjects, obj)
                end
            end
        end
    end
    visualize!(scenes, fobjects;kvs...)
end

function visualize(fig::Figure, objects::Tuple{Any, Vararg{Any}};AxisType=LScene, kwargs...)
    if isa(AxisType, AbstractVector)
        scenes = Any[]
        for _AxisType in AxisType
            lscene = create_axis(_AxisType, fig;kwargs...)
            push!(scenes, lscene)
        end
    else
        lscene = create_axis(AxisType, fig;kwargs...)
        scenes = [lscene for _ in 1:length(objects)]
    end
    visualize!(scenes, objects;kwargs...)
end

function visualize!(scenes::AbstractVector, objects;kwargs...)
    for (lscene,obj) in zip(scenes,objects)
        visualize!(lscene, obj;kwargs...)
    end
end

struct ViewRepresentation <: AbstractRepresentation{Float32, Float64}
    gaze::Vector{Vector{Point3f}}
    position::Vector{Vector{Point2f}}
    timestamp::Vector{Vector{Float64}}
    time_window::Vector{Vector{Float64}}
    event::Vector{Vector{Float64}}
end

function ViewRepresentation(gaze::Vector{Vector{Point3f}}, timestamp::Vector{Vector{Float64}}, time_window, event)
    position = [[Point2f(NaN) for _ in length(gaze[i])] for i in 1:length(gaze)]
    ViewRepresentation(gaze, position, timestamp, time_window, event)
end

get_rep(vp::ViewRepresentation) = vp.position

function ViewRepresentation(spikes::Spiketrain, rp::RippleData, gdata::Union{GazeOnMaze,UnityRaytraceData};fixations_only=true, kwargs...)
    nt = numtrials(gdata)
    gaze = Vector{Vector{Point3f}}(undef, nt)
    pos = Vector{Vector{Point2f}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    timestamp = Vector{Vector{Float64}}(undef, nt)
    time_window = Vector{Vector{Float64}}(undef, nt)

    sp = spikes.timestamps/1000.0 #convert to seconds
    for i in 1:nt
        tg,_gaze,_pos,fixmask,fo = get_trial(gdata,i)
        if length(tg) == 0
            gaze[i] = Point3f[]
            pos[i] = Point2f[]
            events[i] = Float64[]
            timestamp[i] = Float64[]
            time_window[i] = Float64[]
            continue
        end
        tg .= tg .- tg[1]
        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])
        # align to trial start
        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        pos[i] = Vector{Point2f}(undef, nspikes)
        gaze[i] = Vector{Point3f}(undef, nspikes)
        events[i] = Vector{Float64}(undef, nspikes)
        timestamp[i] = Vector{Float64}(undef, nspikes)
        time_window[i] = Vector{Float64}(undef, nspikes)
        js = 1
        for j in 1:nspikes
            k = searchsortedlast(tg,sp_trial[j])
            if 0 < k < size(_gaze,2) && (fixmask[k] || !fixations_only)
                if fo[k] in ["HintImage","CueImage"] 
                    continue
                end
                gaze[i][js] = Point3f(_gaze[:,k])
                pos[i][js] = Point2f(_pos[:,k])
                events[i][js] = sp_trial[j]
                timestamp[i][js] = tg[k]
                time_window[i][js] = tg[k+1]-tg[k]
                js += 1
            end
        end
        pos[i] = pos[i][1:js-1]
        gaze[i] = gaze[i][1:js-1]
        events[i] = events[i][1:js-1]
        timestamp[i] = timestamp[i][1:js-1]
        time_window[i] = time_window[i][1:js-1]
    end
    ViewRepresentation(gaze, pos, timestamp, time_window, events)
end

function ViewRepresentation(gaze_type::Type{T};kwrgas...) where T <: Union{GazeOnMaze, UnityRaytraceData}
    gdata = cd(DPHT.process_level(T)) do
        T()
    end
    rp = cd(DPHT.process_level(RippleData)) do
        RippleData()
    end
    sp = Spiketrain()
    ViewRepresentation(sp,rp,gdata)
end

numtrials(vrp::ViewRepresentation) = length(vrp.position)

function Makie.convert_arguments(::Type{<:AbstractPlot}, vr::ViewRepresentation)
    gazepos = Point3f[]
    for pp in vr.position
        for pq in pp
            push!(gazepos, pq)
        end
    end
    ax3 = S.Axis3(plots=[S.Scatter(gazepos)])
    S.GridLayout(ax3)
end

function create_axis(obj::ViewRepresentation, fig;kwargs...)
    axtype = get_axis_type(ViewRepresentation)
    ax = create_axis(axtype,fig;kwargs...)
end

function visualize!(lscene, vrp::ViewRepresentation;trial::Observable{Trial}=Observable(Trial(1)),kwargs...)
    nt = numtrials(vrp)
    gaze_pos = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return vrp.position[_trial.i]
        else
            return [Point3f(NaN)]
        end
    end
    scatter!(lscene, gaze_pos, color=:green)
end

struct JointRepresentation <: AbstractRepresentation{Float32, Float64}
    events::Vector{Vector{Float64}}
    data::Vector{Matrix{Float32}}
end

get_rep(jp::JointRepresentation) = jp.data

function JointRepresentation(spikes::Spiketrain, rp::RippleData, gdata::Union{GazeOnMaze,UnityRaytraceData}, udata::UnityData)
    sp = spikes.timestamps/1000.0 
    nt = numtrials(gdata)
    nt == numtrials(udata) || error("`gdata` and `udata` should have the same number of trials")
    events = Vector{Vector{Float64}}(undef, nt)
    data = Vector{Matrix{Float32}}(undef, nt)
    for i in 1:nt
        tg,gaze,fixmask = get_trial(gdata,i)
        if isempty(tg)
            events[i] = Float64[]
            data[i] = Matrix{Float32}(undef, 0,0)
            continue
        end
        tg .-= tg[1]
        tu,posx,posy,hd = get_trial(udata,i)
        tu .-= tu[1]

        timestamps = rp.timestamps[i,:]
        
        # find the index of of each spike in this trial
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])

        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        trialdata = zeros(Float32, 5, nspikes)
        trialevents = zeros(Float64, nspikes)
        js = 1
        for j in 1:nspikes
            kg = searchsortedfirst(tg,sp_trial[j])
            ku = searchsortedfirst(tu,sp_trial[j])

            if (0 < kg <= size(gaze,2) && fixmask[kg]) && (0 < ku <= length(posx))
                trialdata[1:3, js] .= gaze[:,kg]
                trialdata[4, js] = posx[ku]
                trialdata[5, js] = posy[ku]
                trialevents[js] = sp_trial[j]
                js += 1
            end
        end
        data[i] = trialdata[:,1:js-1]
        events[i] = trialevents[1:js-1]
    end
    JointRepresentation(events, data)
end

function JointRepresentation(::Type{T};kwargs...) where T <: Union{GazeOnMaze, UnityRaytraceData}
    gdata = cd(DPHT.process_level(T)) do
        T(;kwargs...)
    end
    rp = cd(DPHT.process_level(RippleData)) do
        RippleData()
    end
    udata = cd(DPHT.process_level(UnityData)) do
        UnityData()
    end
    sp = Spiketrain()
    JointRepresentation(sp,rp,gdata,udata)
end


abstract type AbstractViewOccupancy end
struct ViewOccupancy <: AbstractViewOccupancy
    counts::Dict
    bins::Dict
    binidx::Vector{Vector{Tuple{Int64,Int64,Int64,Symbol}}}
    mm::MazeModel
end

struct ViewOccupancyNew{T<:Real} <: AbstractViewOccupancy
    weight::Vector{T}
    mm::SimpleMesh
end

struct ViewAndPlaceOccupancy{T<:Real}
    weight_place::Matrix{T}
    placebin_idx::Vector{Vector{Int64}}
    weight_view::Array{T,3}
    viewbin_idx::Vector{Vector{Int64}}
    mm::SimpleMesh
end

DPHT.level(::Type{<:AbstractViewOccupancy}) = "session"
DPHT.filename(::Type{ViewOccupancy}) = "view_occupancy.mat"
DPHT.filename(::Type{ViewOccupancyNew{T}}) where T <: Real = "view_occupancy_new.jld2"
DPHT.filename(::Type{ViewOccupancyNew}) = "view_occupancy_new.jld2"
DPHT.filename(::Type{ViewAndPlaceOccupancy}) = "view_and_occupancy.jld2"
DPHT.filename(::Type{ViewAndPlaceOccupancy{T}})  where T <: Real = "view_and_occupancy.jld2"
DPHT.level(::Type{ViewAndPlaceOccupancy}) = "session"
DPHT.level(::Type{ViewAndPlaceOccupancy{T}})  where T <: Real = "session"

struct JointOccupancy{T<:Real}
    weight::Dict{CartesianIndex{4},T}
    index::Vector{Vector{CartesianIndex{3}}} # view × place × hd
end

DPHT.filename(::Type{JointOccupancy}) = "joint_occupancy.jld2"
DPHT.filename(::Type{JointOccupancy{T}})  where T <: Real = "joint_occupancy.jld2"
DPHT.level(::Type{JointOccupancy}) = "session"
DPHT.level(::Type{JointOccupancy{T}})  where T <: Real = "session"

function load_jld2(::Type{ViewOccupancyNew})
    fname = DPHT.filename(ViewOccupancyNew)
    fname = replace(fname, ".mat"=>".jld2")
    meta,data = JLD2.load(fname, "meta","data")
    mm = get_maze_mesh()
    ViewOccupancyNew(data["counts"],mm)
end

function save_jld2(voc::ViewOccupancyNew{T};append_tag=true) where T <: Real
    fname = DPHT.filename(ViewOccupancyNew{T})
    fname = replace(fname, ".mat"=>".jld2")
    metadata = Dict{String,Any}() 
    if append_tag
        tag!(metadata, storepatch=true)
    end
    JLD2.save(fname, Dict("data"=>Dict("weight"=>voc.weight), "meta"=>metadata))
end

function save_jld2(vpp::ViewAndPlaceOccupancy{T};append_tag=true) where T <: Real
    fname = DPHT.filename(ViewAndPlaceOccupancy)
    # convert to dictionary
    weight_view_size = size(vpp.weight_view)
    nzidx = findall(vpp.weight_view .!= 0)
    data = Dict("weight_place"=>vpp.weight_place,
                "placebin_idx"=>vpp.placebin_idx,
                "weight_view"=>(idx=nzidx, val=vpp.weight_view[nzidx],vsize=weight_view_size),
                "viewbin_idx"=>vpp.viewbin_idx)
    metadata = Dict{String,Any}() 
    if append_tag
        tag!(metadata, storepatch=true)
    end
    JLD2.save(fname, Dict("data"=>data, "meta"=>metadata))
end

function load_jld2(::Type{ViewAndPlaceOccupancy})
    fname = DPHT.filename(ViewAndPlaceOccupancy)
    meta,data = JLD2.load(fname, "meta", "data")
    weight_view_size = data["weight_view"].vsize
    weight_view_values = data["weight_view"].val
    weight_view_idx = data["weight_view"].idx
    weight_view = zeros(eltype(weight_view_values), weight_view_size...)
    weight_view[weight_view_idx] .= weight_view_values
    mm = get_maze_mesh()
    T = eltype(data["weight_place"])
    ViewAndPlaceOccupancy{T}(data["weight_place"],data["placebin_idx"], weight_view,
                             data["viewbin_idx"], mm)
end

function ViewOccupancy(gdata::Union{GazeOnMaze,UnityRaytraceData}, mm::MazeModel)
    counts,bins,idx = compute_histogram(gdata,mm)
    ViewOccupancy(counts,bins,idx, mm)
end

function ViewOccupancyNew(gdata::Union{GazeOnMaze,UnityRaytraceData}, mm::SimpleMesh;fixations_only=false)
    ng = sum(length.(gdata.gaze) .-1 )
    gaze = zeros(3, ng)
    weight = zeros(ng)
    offset = 0
    for (_gaze, _timestamps,_fix) in zip(gdata.gaze, gdata.timestamps,gdata.fixating)
        Δt = diff(_timestamps)
        nb = length(Δt)
        if fixations_only
            weight[offset+1:offset+nb] = Δt.*_fix[1:end-1]
            gaze[:,offset+1:offset+nb] = _gaze[:,1:end-1] # Skip the last point
        else
            weight[offset+1:offset+nb] = Δt
            gaze[:,offset+1:offset+nb] = _gaze[:,1:end-1] # Skip the last point
        end
        offset += nb
    end
    counts = count_on_manifold(mm, gaze, weight)
    ViewOccupancyNew(counts,mm)
end

function ViewAndPlaceOccupancy(gdata::Union{GazeOnMaze,UnityRaytraceData}, udata::UnityData, mm::SimpleMesh;fixations_only=false, check_dist=false)
    nt = numtrials(gdata)
    nt == numtrials(udata) || error("View and gaze data should have the same number of trials")
    ss = Slice(x=(-12.5, 12.5), y=(-12.5, 12.5), z=(0.0, 0.0))
    m_floor = ss(mm)
    kn = KNearestSearch(mm,1)
    kn_floor = KNearestSearch(m_floor,1)
    weight_place = zeros(nelements(m_floor))
    placebin_idx = Vector{Vector{Int64}}(undef, nt)
    viewbin_idx = Vector{Vector{Int64}}(undef, nt)
    weight_view = zeros(nelements(mm), size(weight_place,1))
    for i in 1:nt
        tg, gaze,_,fixmask,fo = get_trial(gdata,i;trial_start=1)
        if isempty(tg)
            placebin_idx[i] = Int64[]
            viewbin_idx[i] = Int64[]
            continue
        end
        fixated_object = fo 
        tg .-= tg[1]
        tu,posx,posy,hd = get_trial(udata, i;trial_start=1)
        tu .-= tu[1]
        Δtu = diff(tu)
        push!(Δtu, maximum(Δtu))
        _placebin_idx = zeros(Int64, length(tu))
        _viewbin_idx = zeros(Int64, length(tg))
        for (j,(px,py,_tu)) in enumerate(zip(posx, posy, tu))
            idx,dd = searchdists(Meshes.Point(px,py,0.0), kn_floor)
            _idx = first(idx)
            _mm = m_floor[_idx]
            Δ=mean(norm.(_mm.vertices .- centroid(_mm)))
            if (dd[1] <= Δ) || (check_dist == false)
                weight_place[_idx] += Δtu[j] 
                _placebin_idx[j] = _idx
            else
                continue
            end
            # find all gaze coords within this place bin 
            idx0 = searchsortedfirst(tg, _tu)
            idx1 = searchsortedlast(tg, _tu+Δtu[j])
            if _idx == 831
                @debug "trialidx" idx0:idx1
            end
            for k in idx0:idx1-1
                # do not include the hint image
                if fixated_object[k] == "HintImage"
                    continue
                end
                pg = gaze[:,k]
                idxv,dd = searchdists(Meshes.Point(pg...),kn)
                _idxv = first(idxv)
                _mm = mm[_idxv]
                Δ=mean(norm.(_mm.vertices .- centroid(_mm)))
                if (dd[1] <= Δ) || (check_dist == false)
                    #debug
                    if _idxv == 3336 && _idx == 831
                        @debug "indices" tg[k+1]-tg[k]
                    end
                    weight_view[_idxv, _idx] += tg[k+1] -tg[k]
                    _viewbin_idx[k] = _idxv
                else
                    @debug "indices" _idx _idxv dd[1] Δ pg
                end
            end
        end
        placebin_idx[i] = _placebin_idx
        viewbin_idx[i] = _viewbin_idx
    end
    ViewAndPlaceOccupancy(weight_place, placebin_idx, weight_view, viewbin_idx, mm)
end

function ViewAndPlaceOccupancy(gdata::UnityRaytraceData, mm::SimpleMesh;fixations_only=false, trial_start=1,kwargs...)
    nt = numtrials(gdata)
    m_floor = Shadow("xy")(floor_topology3())
    kn = KNearestSearch(mm,1)
    kn_floor = KNearestSearch(m_floor,1)
    weight_place = zeros(Float64,nelements(m_floor),nt)
    placebin_idx = Vector{Vector{Int64}}(undef, nt)
    viewbin_idx = Vector{Vector{Int64}}(undef, nt)
    weight_view = zeros(Float64, nelements(mm), size(weight_place,1),nt)
    for i in 1:nt
        tt, gaze,pos,fixmask,fo = get_trial(gdata,i;trial_start=1)
        if isempty(tt)
            placebin_idx[i] = Int64[]
            viewbin_idx[i] = Int64[]
            continue
        end
        Δt = diff(tt)
        push!(Δt, maximum(Δt))
        _placebin_idx = zeros(Int64, length(tt))
        _viewbin_idx = zeros(Int64, length(tt))
        for (j,(_pos, _gaze, _fo)) in enumerate(zip(eachcol(pos), eachcol(gaze), fo))
            if _fo in ["HintImage","CueImage"]
                continue
            end
            _idx, _idxv = (0,0)
            px,py = _pos[1:2] # don't use the z-coordinate here
            # place bin
            idx,dd = searchdists(Meshes.Point(px,py), kn_floor)
            _idx = first(idx)
            _mm = m_floor[_idx]
            Δ = mean(norm.(_mm.vertices .- centroid(_mm)))
            if dd[1] <= Δ
                weight_place[_idx,i] += Δt[j] 
                _placebin_idx[j] = _idx
            else
                _idx = 0
                continue
            end

            idxv,dd = searchdists(Meshes.Point(_gaze...), kn)
            _idxv = first(idxv)
            _mm = mm[_idxv]
            Δ = mean(norm.(_mm.vertices .- centroid(_mm)))
            if dd[1] <= Δ
                weight_view[_idxv,_idx,i] += Δt[j] 
                _viewbin_idx[j] = _idxv
            else
                continue
            end
        end
        placebin_idx[i] = _placebin_idx
        viewbin_idx[i] = _viewbin_idx
    end
    ViewAndPlaceOccupancy(weight_place, placebin_idx, weight_view, viewbin_idx, mm)
end

function process_kwargs(::Type{JointOccupancy};trial_start=1, nrefinements=(p=3, g=3),kwargs...)
    h = zero(UInt32)
    if trial_start != 1
        h = crc32c(string((:trial_start=>trial_start)),h)
    end
    if nrefinements.p != 3
        h = crc32c(string((:nrefinements_p=>nrefinements.p)),h)
    end
    if nrefinements.g != 3
        h = crc32c(string((:nrefinements_g=>nrefinements.g)),h)
    end
    h
end

function JointOccupancy(gdata::UnityRaytraceData;trial_start=1,nrefinements=(p=3, g=3),kwargs...)
    nt = numtrials(gdata)
    hd_bins = range(0.0, stop=2π, length=24)
    m_floor = Shadow("xy")(floor_topology3(;nrefinements=nrefinements.p))
    kn_floor = KNearestSearch(m_floor,1)
    mm = get_maze_mesh(;nrefinements=nrefinements.g)
    kn = KNearestSearch(mm,1)
    weight = Dict{CartesianIndex{4},Float64}()
    aindex = Vector{Vector{CartesianIndex{3}}}(undef, nt)
    for i in 1:nt
        tt, gaze,pos,fixmask,fo,hd = get_trial(gdata,i;trial_start=trial_start)
        if isempty(tt)
            aindex[i] = CartesianIndex{3}[]
            continue
        end
        Δt = diff(tt)
        push!(Δt, maximum(Δt))
        aindex[i] = [CartesianIndex(0,0,0) for _ in 1:length(tt)]
         for (j,(_pos, _gaze, _fo, _hd)) in enumerate(zip(eachcol(pos), eachcol(gaze), fo,hd))

            if _fo in ["HintImage","CueImage"]
                continue
            end
            # do head direction first
            hidx = argmax(cos.(_hd .- hd_bins))
            _idx, _idxv = (0,0)
            px,py = _pos[1:2] # don't use the z-coordinate here
            # place bin
            idx,dd = searchdists(Meshes.Point(px,py), kn_floor)
            _pidx = first(idx)
            _mm = m_floor[_pidx]
            Δ = mean(norm.(_mm.vertices .- centroid(_mm)))
            if dd[1] > Δ
                continue
            end

            idxv,dd = searchdists(Meshes.Point(_gaze...), kn)
            _idxv = first(idxv)
            _mm = mm[_idxv]
            Δ = mean(norm.(_mm.vertices .- centroid(_mm)))
            if dd[1] > Δ
                continue
            end
            qq = CartesianIndex(_idxv, _pidx, hidx, i)
            weight[qq] = get(weight, qq, 0.0) + Δt[j]
            aindex[i][j] = CartesianIndex(_idxv, _pidx, hidx)
        end
    end
    JointOccupancy(weight,aindex)
end

function JointOccupancy(;redo=false, do_save=true,kwargs...)
    fname = DPHT.filename(JointOccupancy)
    h = process_kwargs(JointOccupancy;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if isfile(fname) && !redo
        jocc = load_jld2(JointOccupancy)
    else 
        unity_gaze_data = UnityRaytraceData(;kwargs...)
        jocc = Hippocampus.JointOccupancy(unity_gaze_data;kwargs...)
        if do_save
            save_jld2(jocc,fname)
        end
    end
    jocc
end

function ViewAndPlaceOccupancy(;do_save=true, redo=false,kwargs...)
    fname = DPHT.filename(ViewAndPlaceOccupancy)
    if isfile(fname) && !redo
        vpp = load_jld2(ViewAndPlaceOccupancy)
    else
        udata = UnityData()
        mm = get_maze_mesh()
        m_floor = floor_topology3()
        unity_gaze_data = UnityRaytraceData()
        vpp = ViewAndPlaceOccupancy(unity_gaze_data, udata, mm,m_floor)
        if do_save
            save_jld2(vpp)
        end
    end
    vpp
end

struct ViewAndPlaceRepresentation <: AbstractRepresentation{Float32,Float64}
    voc::ViewAndPlaceOccupancy{Float64}
    events::Vector{Vector{Float64}}
    placeidx::Vector{Vector{Int64}}
    viewidx::Vector{Vector{Int64}}
end

DPHT.level(::Type{ViewAndPlaceRepresentation}) = "cell"
DPHT.level(::ViewAndPlaceRepresentation) = "cell"

function ViewAndPlaceRepresentation(spikes::Spiketrain, rp::RippleData, udata::UnityData, gdata::UnityRaytraceData, voc::ViewAndPlaceOccupancy{T};fixation_only=false) where T <: Real
    sp = spikes.timestamps/1000.0 
    nt = numtrials(gdata)
    nt == numtrials(udata) || error("`gdata` and `udata` should have the same number of trials")
    events = Vector{Vector{Float64}}(undef, nt)
    viewidx = Vector{Vector{Int64}}(undef, nt)
    placeidx = Vector{Vector{Int64}}(undef, nt)
    for i in 1:nt
        tg,gaze,pos, fixmask,fo = get_trial(gdata,i)
        if isempty(tg)
            events[i] = Float64[]
            viewidx[i] = Int64[]
            placeidx[i] = Int64[]
            continue
        end
        tg .-= tg[1]
        tu,posx,posy,hd = get_trial(udata,i)
        tu .-= tu[1]

        timestamps = rp.timestamps[i,:]
        
        # find the index of of each spike in this trial
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])

        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        trialevents = zeros(Float64, nspikes)
        _viewidx = zeros(Int64, nspikes)
        _placeidx = zeros(Int64, nspikes)
        js = 1
        for j in 1:nspikes
            kg = searchsortedlast(tg,sp_trial[j])
            ku = searchsortedlast(tu,sp_trial[j])
            # check view and place bins
            if (0 < kg <= size(gaze,2) && (fixmask[kg] || !fixation_only)) && (0 < ku <= length(posx))
                #if (voc.placebin_idx[i][ku] != 0) && (voc.viewbin_idx[i][kg] != 0)
                trialevents[js] = sp_trial[j]
                _viewidx[js] = voc.viewbin_idx[i][kg]
                _placeidx[js] = voc.placebin_idx[i][ku]
                js += 1
                #end
            end
        end
        events[i] = trialevents[1:js-1]
        viewidx[i] = _viewidx[1:js-1]
        placeidx[i] = _placeidx[1:js-1]
    end
    ViewAndPlaceRepresentation(voc, events, placeidx, viewidx)
end


struct ViewAndPlaceRepresentationNew <: AbstractRepresentation{Float32,Float64}
    events::Vector{Vector{Float64}}
    placeviewidx::Vector{Vector{Int64}}
end

DPHT.level(::Type{ViewAndPlaceRepresentationNew}) = "cell"
DPHT.level(::ViewAndPlaceRepresentationNew) = "cell"

function ViewAndPlaceRepresentationNew(spikes::Spiketrain, rp::RippleData, gdata::UnityRaytraceData;fixation_only=false)
    sp = spikes.timestamps/1000.0 
    nt = numtrials(gdata)
    events = Vector{Vector{Float64}}(undef, nt)
    placeviewidx = Vector{Vector{Int64}}(undef, nt)
    for i in 1:nt
        tg,gaze,pos, fixmask,fo = get_trial(gdata,i)
        if isempty(tg)
            events[i] = Float64[]
            placeviewidx[i] = Int64[]
            continue
        end
        tg .-= tg[1]
        timestamps = rp.timestamps[i,:]
        
        # find the index of of each spike in this trial
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])

        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        trialevents = zeros(Float64, nspikes)
        _placeviewidx = zeros(Int64, nspikes)
        js = 1
        for j in 1:nspikes
            kg = searchsortedlast(tg,sp_trial[j])
            if (0 < kg <= size(gaze,2) && (fixmask[kg] || !fixation_only))
                #if (voc.placebin_idx[i][ku] != 0) && (voc.viewbin_idx[i][kg] != 0)
                trialevents[js] = sp_trial[j]
                _placeviewidx[js] = kg
                js += 1
                #end
            end
        end
        events[i] = trialevents[1:js-1]
        placeviewidx[i] = _placeviewidx[1:js-1]
    end
    ViewAndPlaceRepresentationNew(events, placeviewidx)
end

function ViewAndPlaceRepresentationNew(;redo=false, do_save=true, kwargs...)
     sp = Spiketrain()
    rp = cd(DPHT.process_level(RippleData)) do
        RippleData()
    end
    unity_gaze_data = cd(DPHT.process_level(UnityRaytraceData)) do
        UnityRaytraceData(;kwargs...)
    end
    vprp = ViewAndPlaceRepresentationNew(sp, rp, unity_gaze_data)
end

function compute_speed(pos::Matrix{T}, timestamp::Array{T}, binidx::Vector{Int64},bmax=maximum(binidx)) where T <: Real
    ds = zero(T)
    dt = zero(T) 
    t0 = timestamp[1] 
    p0 = pos[:,1]
    vv = zeros(bmax)
    current_b = binidx[1] 
    for (p,t,b) in zip(eachcol(pos)[2:end], timestamp[2:end], binidx[2:end])
        if b == 0
            continue
        end
        if b == current_b
            ds += norm(p-p0)
            dt += t-t0
        else
            if current_b > 0
                vv[current_b] += ds/dt
            end
            ds = zero(T)
            dt = zero(T)
            current_b = b
        end
        t0 = t
        p0 .= p
    end
    vv
end

function compute_speed(pos::Vector{Matrix{T}}, timestamp::Vector{Vector{T}}, binidx::Vector{Vector{Int64}},bmax=maximum(maximum.(binidx))) where T <: Real
    vv = zeros(bmax, length(timestamp))
    for i in 1:size(vv,2)
        if isempty(timestamp[i])
            continue
        end
        vv[:,i] = compute_speed(pos[i], timestamp[i], binidx[i], bmax)
    end
    vv
end

function ViewAndPlaceRepresentation(;redo=false,do_save=true)
    sp = Spiketrain()
    rp = cd(DPHT.process_level(RippleData)) do
        RippleData()
    end
    udata = cd(DPHT.process_level(UnityData))  do
        UnityData()
    end
    vpp = cd(DPHT.process_level(ViewAndPlaceOccupancy)) do
        ViewAndPlaceOccupancy()
    end
    unity_gaze_data = cd(DPHT.process_level(UnityRaytraceData)) do
        UnityRaytraceData()
    end
    vprp = ViewAndPlaceRepresentation(sp, rp,udata, unity_gaze_data, vpp)
end

function create_path(posx::AbstractVector{T}, posy::AbstractVector{T}) where T <: Real
    p = [(posx[1], posy[1])]
    for j in 2:length(posx)
        pn = (posx[j], posy[j])
        if pn !== p[end]
            push!(p, pn)
        end
    end
    p 
end

function explore(vvp::ViewAndPlaceOccupancy{<:Real},idxpath::Union{Vector{Int64},Nothing}=nothing;floor_offset=-10.0)
    mm = vvp.mm
    ss = Slice(x=(-12.5, 12.5), y=(-12.5, 12.5), z=(0.0, 0.0))
    m_floor = ss(mm) 
    # this is kind of hacking; find the minimum distance between elements of the floor
    Δ = minimum(norm.(diff(centroid.(m_floor))))
    kn_floor = KNearestSearch(m_floor,1) 
    fig = Figure()
    lscene = LScene(fig[1,1])

    cc1 = cameracontrols(lscene.scene)
    lookat0 = cc1.lookat[]
    eyepos0 = cc1.eyeposition[]
    upvector0 = cc1.upvector[]
    #cc0 = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
    ii = Observable(1)

    _tcolor = fill(NaN, size(vvp.weight_view,1))
    _tcolor .= vvp.weight_view[:,1]
    _alpha = fill(1.0, length(_tcolor))
    _alpha[_tcolor.==0] .= 0.0
    tcolor = Observable(_tcolor)
    alpha = Observable(_alpha)
    pp = coords(centroid(m_floor[ii[]]))

    point = Observable(Point3f(pp.x.val, pp.y.val, pp.z.val))
    fpoint = Observable(Point3f(pp.x.val, pp.y.val, pp.z.val))
    fp = false
    on(ii) do _ii
        if 0 < _ii <= size(vvp.weight_view,2)
            _tcolor  .= vvp.weight_view[:,_ii]
            fill!(_alpha, 1.0)
            _alpha[_tcolor.==0] .= 0.0
            tcolor[] = _tcolor
            alpha[] = _alpha
            pp = coords(centroid(m_floor[ii[]]))
            point[] = Point3f(pp.x.val, pp.y.val, 1.0)
            fpoint[] = Point3f(pp.x.val, pp.y.val, floor_offset)
            if fp
                # also update the camera
                cc1.eyeposition[] = point[]
                cc1.lookat[] = point[] + Makie.Vec(1.0, 0.0, 0.0)
                cc1.upvector[] = Makie.Vec(0.0, 0.0, 1.0)
                update_cam!(lscene.scene, cc1)
            end
        end
    end
    ii[] = 1
    viz!(lscene, vvp.mm;showsegments=true,color=tcolor,alpha=alpha)
    if idxpath !== nothing
        # also plot the path
        pth = centroid.(mm[idxpath])
        pth_points = [Point3f(p.coords.x.val, p.coords.y.val, p.coords.z.val) for p in pth]
        scatter!(lscene, pth_points, color=:black)
    end
    #m_floor2 = Meshes.Translate(0.0, 0.0, floor_offset)(Shadow("xy")(m_floor))

    #viz!(lscene, m_floor2;showsegments=true, alpha=0.0)
    scatter!(lscene, point, color=:red)
    scatter!(lscene,fpoint,color=:red)
    j = 1
    on(events(lscene.scene).keyboardbutton, priority=20) do event
            if ispressed(lscene.scene, Keyboard.c)
                fp = ~fp
                if fp
                    cc1.eyeposition[] = point[]
                    cc1.lookat[] = point[] + Makie.Vec(1.0, 0.0, 0.0)
                    cc1.upvector[] = Makie.Vec(0.0, 0.0, 1.0)
                else
                    cc1.eyeposition[] = eyepos0 
                    cc1.lookat[] = lookat0
                    cc1.upvector[] = upvector0
                end
                update_cam!(lscene.scene, cc1)
            end
            if idxpath === nothing
                μ = centroid(mm[ii[]])
                v = (0.0, 0.0, 0.0).*Unitful.m
                if ispressed(lscene.scene, Keyboard.up)
                    v = (0.0*Unitful.m, Δ, 0.0*Unitful.m)
                elseif ispressed(lscene.scene, Keyboard.down)
                    v = (0.0*Unitful.m, -Δ, 0.0*Unitful.m)
                elseif ispressed(lscene.scene, Keyboard.left)
                    v = (-Δ, 0.0*Unitful.m, 0.0*Unitful.m)
                elseif ispressed(lscene.scene, Keyboard.right)
                    v = (Δ, 0.0*Unitful.m, 0.0*Unitful.m)
                end
                μ1 = Meshes.Translate(v...)(μ)
                _idx,dd = searchdists(μ1, kn_floor)
                if dd[1] <= Δ
                    ii[] = first(_idx)
                end
            else
                if ispressed(lscene.scene, Keyboard.up)
                    j = min(j+1, length(idxpath))
                elseif ispressed(lscene.scene, Keyboard.down)
                    j = max(j-1, 1)
                end
                ii[] = idxpath[j]
            end

        end
    fig
end

function ViewOccupancy(;do_save=true, redo=false)
    fname = DPHT.filename(ViewOccupancy)
    if !redo && isfile(fname)
        vo = DPHT.load(ViewOccupancy)
    else
        mm = MazeModel()
        gdata = UnityRaytraceData()
        vo = ViewOccupancy(gdata,mm)
        if do_save
            DPHT.save(vo)
        end
    end
    vo
end

function ViewOccupancyNew(;do_save=true, redo=false,kwargs...)
    fname = DPHT.filename(ViewOccupancyNew)
    if !redo && isfile(fname)
        vo = DPHT.load(ViewOccupancyNew)
    else
        points,cnx = maze_topology3()
        mm = SimpleMesh(points, connect.(cnx))
        mm2 = refine(refine(refine(mm, QuadRefinement()), QuadRefinement()),QuadRefinement())
        mm2 = SimpleMesh(mm2.vertices, convert(HalfEdgeTopology, mm2.topology))
        gdata = UnityRaytraceData()
        vo = ViewOccupancyNew(gdata,mm2;kwargs...)
        if do_save
            save_jld2(vo)
        end
    end
    vo
end

function Makie.convert_arguments(voc::ViewOccupancy, mm::MazeModel)
    
end

function create_axis(obj::ViewOccupancy, fig;kwargs...)
    axtype = get_axis_type(ViewOccupancy)
    ax = create_axis(axtype,fig;kwargs...)
end

function visualize!(lscene, voc::ViewOccupancy;kwargs...)
    colors = get_maze_colors(voc.mm,voc.counts;kwargs...)
    visualize!(lscene, voc.mm;color=colors,kwargs...)
end


struct ViewMapOld
    counts::Dict
    bins::Dict
    occupancy::Dict
    mm::MazeModel
end

struct ViewMap{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    zbins::AbstractVector{T}
    weight::Array{T,3}
    occupancy::Array{T,3}
end

struct ViewMapNew{T<:Real} <: AbstractMap
    mm::SimpleMesh
    weight::Vector{T}
    occupancy::Vector{T}
end

DPHT.filename(::Type{ViewMapNew}) = "view_map.jld2"
DPHT.filename(::Type{ViewMapNew{T}}) where T <: Real = "view_map.jld2"
DPHT.level(::Type{ViewMapNew}) = "cell"
DPHT.level(::Type{ViewMapNew{T}}) where T <: Real = "cell"

struct SmoothedViewMap{T<:Real}
    mm::SimpleMesh
    weight::Vector{T}
    occupancy::Vector{T}
    unoccupied::Vector{Int64}
    α::T
end

function ViewMapNew(vrp::ViewRepresentation, voc::ViewOccupancyNew{T}) where T <: Real
    gaze = get_gaze(vrp)
    Z = count_on_manifold(voc.mm, gaze)
    ViewMapNew{T}(voc.mm, Z, voc.weight)
end

function process_kwargs(::Type{<:AbstractMap};min_place_duration=0.05, min_view_duration=0.01, min_view_obs=5, min_place_obs=5,kwargs...)
     h = UInt32(0)
    # only store these if they are different from the default
    if min_place_duration != 0.05
        h = crc32c(string(:min_place_duration=>min_place_duration),h)
    end
    if min_view_duration != 0.01
        h = crc32c(string(:min_view_duration=>min_view_duration),h)
    end
    if min_view_obs != 5
        h = crc32c(string(:min_view_obs=>min_view_obs),h)
    end
    if min_place_obs != 5
        h = crc32c(string(:min_place_obs=>min_place_obs),h)
    end
    h
end

function ViewMapNew(vrp::ViewRepresentation, vpoc::ViewAndPlaceOccupancy{T};min_place_duration=0.05, min_view_duration=0.01, min_place_obs=5, min_view_obs=5) where T <: Real
    h = process_kwargs(ViewMapNew;min_place_duration=min_place_duration, min_place_obs=min_place_obs, min_view_obs=min_view_obs)
    gaze = get_gaze(vrp)
    good_place_bins = findall(dropdims(sum(dropdims(sum(vpoc.weight_view,dims=1),dims=1) .> min_place_duration,dims=2),dims=2) .> min_place_obs)
    good_view_bins = findall(dropdims(sum(dropdims(sum(vpoc.weight_view[:,good_place_bins,:],dims=2),dims=2) .> min_view_duration,dims=2),dims=2) .> min_view_obs)
    Z = count_on_manifold(vpoc.mm, gaze)
    # sum over place to get view occupancy
    Zo = dropdims(sum(vpoc.weight_view,dims=(2,3)),dims=(2,3))
    occupancy = zero(Zo)
    occupancy[good_view_bins] .= Zo[good_view_bins]
    Zg = zero(Z)
    Zg[good_view_bins] .= Z[good_view_bins]
    ViewMapNew{T}(vpoc.mm, Zg, occupancy),h
end

function adaptive_smoothing(vm::ViewMapNew{T}, α=T(10000.0)^2;filter_unoccupied=true) where T <: Real
    unoccupied = findall(vm.voc.counts.==0)
    Z,X,Y = adaptive_smoothing(vm.weight, vm.voc.counts, vm.voc.mm, α)
    SmoothedViewMap(vm.voc.mm, X, Y, unoccupied,α)
end

ViewMap(xbins,ybins, zbins, weights) = ViewMap(xbins, ybins, zbins, weights, ones(eltype(weights), size(weights)...))


struct ViewAndPlaceMap{T<:Real} <: AbstractMap
    mm::SimpleMesh
    weight::Matrix{T}
    occupancy::Matrix{T}
end

DPHT.filename(::Type{ViewAndPlaceMap}) = "view_and_place_map.jld2"
DPHT.filename(::ViewAndPlaceMap{T}) where T <: Real = "view_and_place_map.jld2"

DPHT.level(::Type{ViewAndPlaceMap}) = "cell"
DPHT.level(::ViewAndPlaceMap{T}) where T <: Real = "cell"


function ViewAndPlaceMap(vprp::ViewAndPlaceRepresentation, vpp::ViewAndPlaceOccupancy{T};smoothing_params::NamedTuple=(method=:gaussian, σ=5)) where T <: Real
    P = zeros(size(vpp.weight_view)...)
    for i in 1:length(vprp.events)
        for j in 1:length(vprp.events[i])
            vidx = vprp.viewidx[i][j]
            pidx = vprp.placeidx[i][j]
            if (pidx > 0) && (vidx > 0)
                P[vprp.viewidx[i][j], vprp.placeidx[i][j]] += 1.0
            end
        end
    end
    occupancy = zeros(T, size(vpoc.weight_view)[1:2]...)
    occupancy[good_view_bins, good_place_bins] .= dropdims(sum(vpoc.weight_view[good_view_bins, good_place_bins,:],dims=3),dims=3)
    ViewAndPlaceMap{T}(vpoc.mm, P, occpancy)
end


function ViewAndPlaceMap(vrp::ViewRepresentation, vpoc::ViewAndPlaceOccupancy{T};min_view_duration=0.01, min_place_duration=0.05, min_view_obs=5, min_place_obs=5) where T <: Real
    good_place_bins = findall(dropdims(sum(dropdims(sum(vpoc.weight_view,dims=1),dims=1) .> min_place_duration,dims=2),dims=2) .> min_place_obs)
    good_view_bins = findall(dropdims(sum(dropdims(sum(vpoc.weight_view[:,good_place_bins,:],dims=2),dims=2) .> min_view_duration,dims=2),dims=2) .> min_view_obs)
    mm = vpoc.mm
    m_floor = Shadow("xy")(floor_topology3())
    #project onto 2D
    m_floor = Shadow("xy")(m_floor)
    pos = get_positions(vrp)
    gaze = get_gaze(vrp)
    kidx_p = mapto(m_floor, Tuple.(eachcol(pos)))
    kidx_g = mapto(mm, Tuple.(eachcol(gaze)))
    P = zeros(T, nelements(mm), nelements(m_floor))
    P2 = zeros(T, nelements(mm), nelements(m_floor))
    for (kg,kp) in zip(kidx_g, kidx_p)
        if !isempty(kg) && !isempty(kp)
            P[kg,kp] .+= one(T)
        end
    end
    P2[good_view_bins, good_place_bins] .= P[good_view_bins, good_place_bins]
    occupancy = zeros(T, size(vpoc.weight_view)[1:2]...)
    occupancy[good_view_bins, good_place_bins] .= dropdims(sum(vpoc.weight_view[good_view_bins, good_place_bins,:],dims=3),dims=3)
    ViewAndPlaceMap{T}(vpoc.mm, P, occupancy)
end

function ViewAndPlaceMap(;redo=false, do_save=true,kwargs...)
    fname = DPHT.filename(ViewAndPlaceMap)
    h = process_kwargs(ViewAndPlaceMap;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2","_$(hs).jld2")
    end
    if isfile(fname) && !redo
        vmpv = load_jld2(ViewAndPlaceMap,fname)
    else
        vrp = cd(DPHT.process_level(ViewAndPlaceRepresentation)) do
            ViewRepresentation(UnityRaytraceData;kwargs...)
        end
        vpoc = cd(DPHT.process_level(ViewAndPlaceOccupancy)) do
            ViewAndPlaceOccupancy(;kwargs...)
        end
        vmpv = ViewAndPlaceMap(vrp, vpoc;kwargs...)
        if do_save
            save_jld2(vmpv,fname)
        end
    end
    vmpv
end

function explore(vmvp::ViewAndPlaceMap{T},ii) where T <: Real
    fig = Figure()
    lg = GridLayout(fig[1,1])
    explore!(lg, vmvp,ii)
    fig
end

function explore!(lg, vmvp::ViewAndPlaceMap{T},ii::Integer) where T <: Real
    tcolor = vmvp.weight_smooth[:,ii]
    tcolor[(!isfinite).(tcolor)] .= 0.0
    alpha = get_alpha(tcolor)
    explore!(lg, vmvp.voc.mm;color=tcolor,alpha=alpha, showsegments=true, floor_offset=-10, ceiling_offset=10, colormap=:jet, label="Firing rate")
end

function save_jld2(vmpv::ViewAndPlaceMap{T};append_tag=true) where T <: Real
    fname = DPHT.filename(ViewAndPlaceMap)
    data = Dict("weight"=>vmpv.weight, "weight_smooth"=>vmpv.weight_smooth, "smoothing_params"=>vmpv.smoothing_params)
    metadata = Dict{String,Any}() 
    if append_tag
        tag!(metadata, storepatch=true)
    end
    JLD2.save(fname, Dict("data"=>data, "meta"=>metadata))
end

function load_jld2(::Type{ViewAndPlaceMap})
    fname = DPHT.filename(ViewAndPlaceMap)
    meta,data = JLD2.load(fname, "meta","data")
    vpp = cd(DPHT.process_level(ViewAndPlaceOccupancy)) do
        ViewAndPlaceOccupancy()
    end
    ViewAndPlaceMap(vpp, data["weight"], data["weight_smooth"], data["smoothing_params"])
end

function ViewMap(vrp::ViewRepresentation, xbins::AbstractVector{T}, ybins::AbstractVector{T}, zbins::AbstractVector{T}) where T <: Real
    view_count = fill(0.0, length(xbins)-1, length(ybins)-1, length(zbins)-1)
    nt = numtrials(vrp)
    for i in 1:nt
        position = vrp.position[i]
        xpos = [pos[1] for pos in position] 
        ypos = [pos[2] for pos in position] 
        zpos = [pos[3] for pos in position] 
        h = fit(Histogram, (xpos,ypos,zpos), (xbins, ybins,zbins))
        view_count .+= h.weights
    end
    ViewMap(xbins,ybins, zbins, view_count)
end

function ViewMap(vrp::ViewRepresentation, voc::ViewOccupancy)
    mm = voc.mm
    bins = get_bins(mm)
    # convert to matrix
    gaze_pos = Vector{Matrix{Float64}}(undef, length(vrp.position))
    for (ii,pos) in enumerate(vrp.position)
        gaze_pos[ii] = fill(0.0, 3, length(pos))
        for (jj,p) in enumerate(pos)
            gaze_pos[ii][:,jj] .= p
        end
    end
    counts, idx = compute_histogram(gaze_pos,bins)
    # normalize using the occupancy map
    for (k,v) in voc.counts
        for (ii,vv) in enumerate(v)
            counts[k][ii] ./= vv
        end
    end

    ViewMap(counts, bins, voc.counts, mm), idx
end

function ViewMap(;kwargs...)
    mm = MazeModel()
    voc = cd(DPHT.process_level(ViewOccupancy)) do
        ViewOccupancy()
    end
    vrp = ViewRepresentation()
    ViewMap(vrp, mm,voc)
end

function ViewMapNew(gaze_type::Type{T};redo=false, do_save=true, kwargs...) where T <: Union{GazeOnMaze, UnityRaytraceData}
    fname = DPHT.filename(ViewMapNew)
    h = process_kwargs(ViewMapNew;kwargs...)
    if h > 0
        hs = string(h, base=16)
        fname = replace(fname, ".jld2"=>"_$(hs).jld2")
    end
    if !redo && isfile(fname)
        vm = load_jld2(ViewMapNew)
    else
        vpoc = cd(DPHT.process_level(ViewAndPlaceOccupancy)) do
            ViewAndPlaceOccupancy(;kwargs...)
        end
        vrp = ViewRepresentation(gaze_type;kwargs...)
        # TODO: Save this
        vm, h = ViewMapNew(vrp, vpoc)
        if do_save
            save_jld2(vm, fname)
        end
    end
    vm
end

function smooth(vm::ViewMap;kwargs...)
    D, points, pidx,ll = compute_distance_matrix(vm.mm)
    Z = smooth(vm, D, pidx)
    Z, D, points, pidx, ll
end

function smooth(vm::ViewMap,D,pidx;kwargs...)
    smooth(vm.counts,D, pidx;kwargs...) 
end

function create_axis(obj::ViewMap, fig;kwargs...)
    axtype = get_axis_type(ViewMap) 
     ax = create_axis(axtype,fig;kwargs...)
 end

function visualize!(lscene, vm::ViewMap;kernel=nothing, normalize=true, kwargs...)
    # normalize each component in vm.counts
    if normalize
        ncounts = typeof(vm.counts)()
        for k in keys(vm.counts)
            cc = vm.counts[k]
            ncounts[k] = Vector{Array{Float64,3}}(undef, length(cc))
            for ii in 1:length(cc)
                ncounts[k][ii] = cc[ii]./vm.occupancy[k][ii]
            end
        end
    else
        ncounts = vm.counts
    end
    colors = get_maze_colors(vm.mm,ncounts;kernel=kernel,kwargs...)
    visualize!(lscene, vm.mm;color=colors,kwargs...)
end
