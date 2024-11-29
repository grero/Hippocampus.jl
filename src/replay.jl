using MAT
using CSV
using DataFrames
using Makie
using ProgressMeter

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

function DPHT.save(gdata::GazeOnMaze)
    fname = DPHT.filename(typeof(gdata))
    qdata = Dict{String,Any}()
    metadata = Dict{String,Any}() 
    tag!(metadata, storepatch=true)
    for k in fieldnames(GazeOnMaze)
        v = getfield(gdata, k)
        qdata[string(k)] = v
    end
    qdata["meta"] = metadata
    MAT.matwrite(fname, qdata)
end

function DPHT.load(::Type{GazeOnMaze})
    fname = DPHT.filename(GazeOnMaze)
    qdata = MAT.matread(fname)
    metadata = qdata["meta"]
    args = Any[]
    for k in fieldnames(GazeOnMaze)
        push!(args, qdata[string(k)])
    end
    GazeOnMaze(args...)
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
    gaze::Vector{Matrix{Float64}}
    position::Vector{Matrix{Float64}}
    head_direction::Vector{Vector{Float64}}
    timestamps::Vector{Vector{Float64}}
end

function UnityRaytraceData()
    edata = cd(DPHT.process_level(EyelinkData)) do
        EyelinkData()
    end
    fname = "unityfile_eyelink.csv"
    if !ispath(fname)
        error("No raytracing data found")
    end
    unity_eyelinkfile = CSV.File(fname, header=0)
    n = length(unity_eyelinkfile)
    fixated_points = fill(NaN, 3, n)
    position = fill(0.0, 3, n)
    direction = fill(0.0, n)
    timestamps = zeros(UInt64,n)
    i = 1
    for row in unity_eyelinkfile
        # TODO: Grab more data here
        px,py,pz,a = (row[6],row[7],row[8],row[9])
        θ = π*a/180.0
        position[:,i] = (px,pz,py)
        direction[i] = θ
        timestamps[i] = row[2]
        gx,gy,gz = (row[10],row[11],row[12])
        if (gx !== missing) && (gy !== missing) && (gz !== missing)
            # gx,gy,gz is relative to player ?

            fixated_points[:,i] .= (gx,gz,gy) # unity has the z-axis into the scene
            i += 1
        end
    end

    # break up into trials using edata
    nt = numtrials(edata)
    trial_fixations = Vector{Matrix{Float64}}(undef,nt)
    trial_position = Vector{Matrix{Float64}}(undef, nt)
    trial_head_direction = Vector{Vector{Float64}}(undef, nt)
    trial_times = Vector{Vector{Float64}}(undef,nt)
    for i in 1:nt
        te,_,_, = get_trial(edata, i)
        idx0 = searchsortedfirst(timestamps, te[1])
        idx1 = searchsortedlast(timestamps,te[end])
        trial_fixations[i] = fixated_points[:,idx0:idx1]
        trial_position[i] = position[:,idx0:idx1]
        trial_head_direction[i] = diration[idx0:idx1]
        trial_times[i] = timestamps[idx0:idx1]/1000.0 # convert to seconds
    end
    UnityRaytraceData(trial_fixations, trial_position, trial_head_diraction, trial_times)
end

function visualize!(lscene, unitygaze::UnityRaytraceData;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0))
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
        tu,tg,tp,td 
    end

    current_pos = Observable(ugdata_trial[][2][1:1])
    current_path = Obserable(ugdata_trial[][3])
    current_dir = Observable(ugdata_trial[][4][1:1])
    current_arrow = lift(current_dir) do θ
        [Point3f(cos(θ), sin(θ), 0.0)]
    end
    current_j = 1

    onany(ugdata_trial, current_time) do _ugt, _ct
        _tu, _points = (_ugt[1], _ugt[2])
        j = searchsortedfirst(_tu, _ct) 
        current_pos[] = _points[current_j:j]
        current_j = j
    end
    lines!(lscene, current_path)
    scatter!(lscene,current_pos, color=:black)
    arrows!(lscene, current_pos, current_arrow, color=:black)
end

"""
Wrapper type to visualize the raytracing
"""
struct RaytraceViewer
    udata::UnityData
    gazemaze::GazeOnMaze
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

"""
Type to indicate that we want to replay the maze as the subject experienced it.
"""
struct MazeReplayer
    mm::MazeModel
    udata::UnityData
end

function visualize!(lscene::LScene, mp::MazeReplayer;current_time::Observable{Float64}=Observable(0.0), trial::Observable{Trial}=Observable(Trial(1)), kwargs...)

    visualize!(lscene,mp.mm;show_ceiling=true)
    cc = cameracontrols(lscene.scene)
    cc.fov[] = 60.0
    cc.near[] = 0.3
    cc.far[] = 1000.0
    cc.settings.clipping_mode[] = :adaptive # :static
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
            cc.lookat[] = Point3f(sin(θ), cos(θ), 0.0) + pos
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

function create_axis(obj::ViewRepresentation, fig;kwargs...)
    axtype = get_axis_type(ViewRepresentation)
    ax = create_axis(axtype,fig;kwargs...)
end

function create_axis(obj::ViewOccupancy, fig;kwargs...)
    axtype = get_axis_type(ViewOccupancy)
    ax = create_axis(axtype,fig;kwargs...)
end

function create_axis(obj::ViewMap, fig;kwargs...)
   axtype = get_axis_type(ViewMap) 
    ax = create_axis(axtype,fig;kwargs...)
end

 get_axis_type(::Type{T}) where T <: Any = LScene
 get_axis_type(::Type{EyelinkData}) = Axis

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

struct ViewRepresentation
    position::Vector{Vector{Point3f}}
    event::Vector{Vector{Float64}}
end

function ViewRepresentation(spikes::Spiketrain, rp::RippleData, gdata::GazeOnMaze)
    nt = numtrials(gdata)
    position = Vector{Vector{Point3f}}(undef, nt)
    events = Vector{Vector{Float64}}(undef, nt)
    sp = spikes.timestamps/1000.0 #convert to seconds
    for i in 1:nt
        tg,gaze,fixmask = get_trial(gdata,i)

        timestamps = rp.timestamps[i,:]
        idx0 = searchsortedfirst(sp, timestamps[1])
        idx1 = searchsortedlast(sp, timestamps[3])
        # align to trial start
        sp_trial = sp[idx0:idx1] .- timestamps[1]
        nspikes = idx1-idx0+1
        position[i] = Vector{Point3f}(undef, nspikes)
        events[i] = Vector{Float64}(undef, nspikes)
        js = 1
        for j in 1:nspikes
            k = searchsortedfirst(tg,sp_trial[j])
            if 0 < k <= size(gaze,2) && fixmask[k]
                position[i][js] = Point3f(gaze[:,k])
                events[i][js] = sp_trial[j]
                js += 1
            end
        end
        position[i] = position[i][1:js-1]
        events[i] = events[i][1:js-1]
    end
    ViewRepresentation(position, events)
end

function ViewRepresentation(;kwrgas...)
    gdata = cd(DPHT.process_level(GazeOnMaze)) do
        GazeOnMaze()
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

struct ViewOccupancy
    counts::Dict
    bins::Dict
    mm::MazeModel
end

DPHT.level(::Type{ViewOccupancy}) = "session"

function ViewOccupancy(gdata::GazeOnMaze, mm::MazeModel)
    counts,bins = compute_histogram(gdata,mm)
    ViewOccupancy(counts,bins, mm)
end

function ViewOccupancy()
    mm = MazeModel()
    gdata = GazeOnMaze()
    ViewOccupancy(gdata,mm)
end

function Makie.convert_arguments(voc::ViewOccupancy, mm::MazeModel)
end

function visualize!(lscene, voc::ViewOccupancy;kwargs...)
    colors = get_maze_colors(voc.mm,voc.counts)
    visualize!(lscene, voc.mm;color=colors,kwargs...)
end


struct ViewMap
    counts::Dict
    bins::Dict
    occupancy::Dict
    mm::MazeModel
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
            gaze_pos[ii][:,jj] = p
        end
    end
    counts = Dict{Symbol,Vector{Array{Float64,3}}}()
    for k in keys(bins)
        counts[k] = compute_histogram(gaze_pos,bins[k])
    end
    ViewMap(counts, bins, voc.counts, mm)
end

function ViewMap(;kwargs...)
    mm = MazeModel()
    voc = cd(DPHT.process_level(ViewOccupancy)) do
        ViewOccupancy()
    end
    vrp = ViewRepresentation()
    ViewMap(vrp, mm,voc)
end

function visualize!(lscene, vm::ViewMap;normalize=true, kwargs...)
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
    colors = get_maze_colors(vm.mm,ncounts)
    visualize!(lscene, vm.mm;color=colors,kwargs...)
end
