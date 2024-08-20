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


function compute_histogram(gdata::GazeOnMaze;fixations_only=true)
    bins,normals = create_maze(;Δz=0.01)
    all_bins = [bins.pillars[:];bins.walls;bins.ceiling;bins.floor]
    if fixations_only 
        gaze = Vector{Matrix{Float64}}(undef, length(gdata.gaze))
        for i in eachindex(gaze)
            gaze[i] = gdata.gaze[i][:,gdata.fixation[i]]
        end
    else
        gaze = gdata.gaze
    end
    counts = compute_histogram(gaze,all_bins)
    counts,bins,normals
end

function show_maze(bins,normals,counts::Union{Vector{Array{T,3}},Nothing}=nothing;explore=false, replay=false, interactive=false, gdata::Union{Nothing, GazeOnMaze}=nothing, udata::Union{Nothing, UnityData}=nothing, trial::Int64=1) where T <: Real
    fig = Figure()
    #ax = Axis3(fig[1,1],aspect=:data)
    lscene = LScene(fig[1,1], show_axis=false)
    for (i,(bin,n)) in enumerate(zip(bins,normals))
        m = CartesianGrid(first.(bin), last.(bin);dims=length.(bin))
        if counts !== nothing
            # we want to color only the inside
            c = counts[i]
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
            _color = RGB(0.8, 0.8, 0.8) 
        end
        viz!(lscene, m, color=_color,colormap=:Blues)
    end

    lookat = Point3f(1.0, 0.0, 2.5)
    if replay
        # replay experiment with the supplied position
        #cc = Makie.Camera3D(lscene.scene, projectiontype = Makie.Perspective, rotation_center=:eyeposition, center=false)
        cc = cameracontrols(lscene.scene)
        cc.fov[] = 39.6
        tg,gaze,fixmask = (gdata.time[trial], gdata.gaze[trial],gdata.fixation[trial])
        tg .-= tg[1]

        tp,posx,posy,head_direction = get_trial(udata,trial)
        tp .-= tp[1]
        position = permutedims([posx posy])
        gaze_pos = Observable([Point3f(NaN)])
        current_j = 1
        ii = Observable(1)

        on(ii) do i
            _tp = tp[i]
            # grab the points 
            j = searchsortedfirst(tg, _tp)
            pos = Point3f(position[1,i], position[2,i], 2.5)
            θ = π*head_direction[i]/180

            # only use fixation points
            _fixmask = fixmask[current_j:j]
            #gaze_pos[] = [Point3f(dropdims(mean(gaze[:,current_j+1:j][:,_fixmask],dims=2),dims=2))]
            gaze_pos[] = Point3f.(eachcol(gaze[:,current_j:j][:,_fixmask]))
            current_j = j
            cc.lookat[] = Point3f(cos(θ), sin(θ), 0.0) + pos
            cc.eyeposition[] = pos
            update_cam!(lscene.scene, cc)
        end
        scatter!(lscene, gaze_pos, color=:red)

        on(events(lscene.scene).scroll, priority=20) do (dx,dy)
            i_new = round(Int64,ii[] + 5*dx)
            if 0 < i_new <= size(position,2)
                ii[] = i_new
            end
        end

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
            tg,gaze,fixmask = (gdata.time[trial], gdata.gaze[trial],gdata.fixation[trial])
            scatter!(lscene, gaze[1,fixmask], gaze[2,fixmask], gaze[3,fixmask],color=:red)
        end
        if udata !== nothing
            tp,posx,posy,dir = get_trial(udata,trial)
            position = permutedims([posx posy])
            lines!(lscene, position[1,:], position[2,:], fill(0.0, size(position,2)),color=:black)
        end
        # show current position
        ii = Observable(1)
        current_pos = lift(ii) do i
            pos = Point3f(position[1,i],position[2,i],0.0)
            [pos]
        end

        on(events(lscene.scene).keyboardbutton, priority=20) do event
            if ispressed(lscene.scene, Keyboard.up)
                ii[] = min(ii[]+1,size(position,2))
            elseif ispressed(lscene.scene, Keyboard.down)
                ii[] = max(ii[]-1,1)
            end
        end
        scatter!(lscene, current_pos, color=:blue)
    end
    fig
end