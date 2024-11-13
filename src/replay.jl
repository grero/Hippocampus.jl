using MAT

struct GazeOnMaze
    time::Vector{Vector{Float64}}
    gaze::Vector{Matrix{Float64}}
    fixation::Vector{Vector{Bool}}
    triggers::Matrix{Int64}
    timestamps::Matrix{Float64}
    header::Dict
end

DPHT.filename(::Type{GazeOnMaze}) = "maze_raytrace.mat"


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
        edata = EyelinkData()
        udata = UnityData()
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

function Makie.convert_arguments(::Type{<:AbstractPlot}, gdata::GazeOnMaze)
    # 3D scatter plot of all positions
    ax3 = S.Axis3(plots=[S.Scatter(gaze[1,:], gaze[2,:], gaze[3,:]) for gaze in gdata.gaze])
    # TODO: Use the "exploded view here", that is indicate the pillars, as well as the floor
    # and ceiling
    S.GridLayout(ax3)
end

function visualize!(lscene, gdata::GazeOnMaze;trial::Observable{Trial}=Observable(Trial(1)), current_time::Observable{Float64}=Observable(0.0))
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
        fixmask = _gdt[3]
        j = searchsortedfirst(tg, _ct) 
        if 0 < j <= length(fixmask)
            _fixmask = fixmask[current_j:j]
            _current_j = current_j
            current_j = j
            current_gaze[] = Point3f.(eachcol(gaze[:,_current_j:j][:,_fixmask]))
        else
            current_gaze[] = [Point3f(NaN)]
        end
    end

    scatter!(lscene, current_gaze, color=RGB(0.8, 0.8, 0.8))
end


function compute_histogram(gdata::GazeOnMaze,mm::MazeModel;fixations_only=true)
    bins = get_bins(mm)
    if fixations_only 
        gaze = Vector{Matrix{Float64}}(undef, length(gdata.gaze))
        for i in eachindex(gaze)
            gaze[i] = gdata.gaze[i][:,gdata.fixation[i]]
        end
    else
        gaze = gdata.gaze
    end
    counts = Dict{Symbol,Vector{Array{Float64,3}}}()
    for k in keys(bins)
        counts[k] = compute_histogram(gaze,bins[k])
    end
    counts,bins,normals
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

function visualize(objects...;kwargs...)
    fig = Figure()
    lscene = LScene(fig[1,1], show_axis=false)
    # attach events
    current_time = Observable(0.0)
    on(events(lscene.scene).scroll, priority=20) do (dx,dy)
        current_time[] = current_time[] + dx
    end
    current_trial = Observable(Trial(1))
    on(events(lscene.scene).keyboardbutton, priority=20) do event
        has_changed = false
        if ispressed(lscene.scene, Keyboard.up)
            current_trial[] = Trial(current_trial[].i+1)
            # TODO: Is there a meaninful way to check whether we have reached the end here?
            has_changed = true
        elseif ispressed(lscene.scene, Keyboard.down)
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
    for obj in objects
        visualize!(lscene,obj;current_time=current_time, trial=current_trial, kwargs...)
    end
    fig
end