using MAT
using CSV
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

get_trial(gdata, i) = (gdata.time[i], gdata.gaze[i], gdata.fixation[i])

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
function raytrace(x, y, pos,direction, focal_length;camera_height=2.5)
    # find the angle of the point
    θ = atan(x, focal_length)
    ϕ = atan(y, focal_length)
    θ += direction
    # now cast along the line θ

    dl = 0.01
    xp,yp,zp = (pos[1], pos[2],camera_height)
    sθ = sin(θ)
    cθ = cos(θ)
    sϕ = sin(ϕ)
    while true
        dx,dy,dz = (dl*cθ, dl*sθ,dl*sϕ)
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
Return the 3D eye position on objects in the maze
"""
function GazeOnMaze(edata::EyelinkData, udata::UnityData)
    # we an only align edata and udata using events, so we need to operate on trials
    # the most compact representation
    nt = size(edata.triggers,1)
    gaze = Vector{Matrix{Float64}}(undef, nt)
    gtime = Vector{Vector{Float64}}(undef,nt)
    fixation = Vector{Vector{Bool}}(undef, nt)
    screen_width = edata.header["gaze_coords"][3] - edata.header["gaze_coords"][1]
    screen_height = edata.header["gaze_coords"][4] - edata.header["gaze_coords"][2]
    gx0, gy0, gxm, gym = edata.header["gaze_coords"]
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
                x = scale_to_camera(gx[j]-gx0, 36.0, screen_width)
                y = scale_to_camera(gy[j]-gy0, 24.0, screen_height)
                # use negative gx and gy here since the camera flips these
                _gaze[:,j] .= raytrace(-x,-y,[posx[tidx],posy[tidx]],π*dir[tidx]/180,50.0)
            end
        end
        gaze[i] = _gaze
        next!(prog)
    end
    GazeOnMaze(gtime,gaze,fixation, edata.triggers, edata.timestamps, Dict("focal_length"=>50.0,"camera_height"=>2.5))
end

function MakieCore.convert_arguments(::Type{<:AbstractPlot}, gdata::GazeOnMaze)
    # 3D scatter plot of all positions
    ax3 = S.Axis3(plots=[S.Scatter(gaze[1,:], gaze[2,:], gaze[3,:]) for gaze in gdata.gaze])
    # TODO: Use the "exploded view here", that is indicate the pillars, as well as the floor
    # and ceiling
    S.GridLayout(ax3)
end

function visualize!(lscene, gdata::GazeOnMaze;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0),fixation_only=true, kwargs...)
    nt = numtrials(gdata)
    current_gaze = Observable([Point3f(NaN)])
    gdata_trial = lift(trial) do _trial
        if 0 < _trial.i <= nt
            return gdata.time[_trial.i], gdata.gaze[_trial.i],gdata.fixation[_trial.i]
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
            _fixmask = fixmask[_current_j:j]
            current_gaze[] = Point3f.(eachcol(gaze[:,_current_j:j][:,_fixmask]))
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

function visualize!(lscene::LScene, mp::MazeReplayer;current_time::Observable{Float64}=Observable(0.0), trial::Observable{Trial}=Observarble(Trial(1)), kwargs...)

    visualize!(lscene,mp.mm)
    cc = cameracontrols(lscene.scene)
    cc.fov[] = 39.6

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
            pos = Point3f(px[j],py[j], 2.5)
            θ = π*dir[j]/180.0 # convert to radians
            cc.lookat[] = Point3f(cos(θ), sin(θ), 0.0) + pos
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

function show_maze!(lscene, bins,counts::Union{Dict{Symbol,Vector{Array{T,3}}},Nothing}=nothing,normals::Union{Dict{Symbol,Vector{Vector{Float64}}},Nothing}=nothing;explore=false, replay=false, interactive=false, gdata::Union{Nothing, GazeOnMaze}=nothing, udata::Union{Nothing, UnityData}=nothing, trial::Observable{Int64}=Observable(1),offsets::Union{Nothing, Dict{Symbol, Vector{Vector{Float64}}}}=nothing,show_ceiling=true,posters=nothing) where T <: Real
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
    
    position = lift(udata_trial) do _udata
        [Point3f(pos[1], pos[2], 0.5) for pos in eachcol(_udata[2])]
    end

    lookat = Point3f(1.0, 0.0, 2.5)
    ii = Observable(1)
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
        cc.fov[] = 39.6
        if gdata !== nothing
            tg,gaze,fixmask = (gdata.time[trial[]], gdata.gaze[trial[]],gdata.fixation[trial[]])
            tg .-= tg[1]
        end

        #tp,posx,posy,head_direction = get_trial(udata,trial[])
        #tp .-= tp[1]
        #position = permutedims([posx posy])
        gaze_pos = Observable([Point3f(NaN)])
        current_j = 1

        on(ii) do i
            # grab the points 
            tp, position, head_direction = udata_trial[]
            _tp = tp[i]
            pos = Point3f(position[1,i], position[2,i], 2.5)
            θ = π*head_direction[i]/180

            # only use fixation points
            #gaze_pos[] = [Point3f(dropdims(mean(gaze[:,current_j+1:j][:,_fixmask],dims=2),dims=2))]
            if gdata !== nothing
                j = searchsortedfirst(tg, _tp)
                gaze_pos[] = Point3f.(eachcol(gaze[:,current_j:j][:,_fixmask]))
                current_j = j
                _fixmask = fixmask[current_j:j]
            end
            cc.lookat[] = Point3f(cos(θ), sin(θ), 0.0) + pos
            cc.eyeposition[] = pos
            update_cam!(lscene.scene, cc)
        end
        scatter!(lscene, gaze_pos, color=:red)

        if !interactive
            @async for j in 1:length(tp)
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
        if gdata !== nothing
            tg,gaze,fixmask = (gdata.time[trial[]], gdata.gaze[trial[]],gdata.fixation[trial[]])
            scatter!(lscene, gaze[1,fixmask], gaze[2,fixmask], gaze[3,fixmask],color=:red)
        end
        # show current position
        ii = Observable(1)
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
        arrows!(lscene, current_pos, current_dir,color=:blue)
    end
end


# TODO: It would be more elegant to make use of Makie recipe here

# I feel like this is duplicating functionality that must be in Makie somwhere...
function create_axis(AxisType,fig;kwargs...)
    if AxisType <: LScene
        axis_args = (show_axis=get(kwargs, :show_axis,false),)
    else
        axis_args = (backgroundcolor=get(kwargs, :backgroundcolor, :white),)
    end
    lscene = AxisType(fig[1,1];axis_args...)
end

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
        if ispressed(fig.scene, Keyboard.up)
            current_trial[] = Trial(current_trial[].i+1)
            # TODO: Is there a meaninful way to check whether we have reached the end here?
            has_changed = true
        elseif ispressed(fig.scene, Keyboard.down)
            nc  = current_trial[].i+1
            if nc > 0
                current_trial[] = Trial(nc)
                has_changed = true
            end
        end
        if has_changed
            current_time[] = 0.0
        end
    end
    visualize(fig, objects;current_time=current_time, trial=current_trial, kwargs...)
    current_trial[] = Trial(1)
    current_time[] = 0.0
    fig
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

function visualize(fig::Figure, objects::Vector{Any};AxisType=[LScene for _ in length(objects)],kwargs...)
    scenes = Any[]
    lg = GridLayout(fig[1,1])
    for (i,(obj,_AxisType)) in enumerate(zip(objects,AxisType))
        if isa(obj, Tuple)
            n = length(obj)
        else
            n = 1
        end
        if isa(_AxisType, NTuple{n,Any})
            lscene = [create_axis(_AT, lg[1,i];kwargs...) for _AT in _AxisType]
        else
            lscene = create_axis(_AxisType, lg[1,i];kwargs...)
        end
        push!(scenes, lscene)
    end
    for (objs,scene) in zip(objects, scenes)
        @show typeof(objs) typeof(scene)
        visualize!(scene, objs;kwargs...)
    end
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

numtrials(vrp::ViewRepresentation) = length(vrp.position)

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

struct ViewMap{T<:Real}
    xbins::AbstractVector{T}
    ybins::AbstractVector{T}
    zbins::AbstractVector{T}
    weights::Array{T,3}
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