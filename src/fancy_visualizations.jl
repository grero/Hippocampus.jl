using Makie
using Meshes
using FileIO

function visualize_trial(qdata::Hippocampus.UnityRaytraceData, udata::Hippocampus.UnityData,vpvrp::Hippocampus.ViewAndPlaceRepresentationNew;start_trial=1, end_trial=nothing,fname::Union{String, Nothing}="test.mp4")
    nt = Hippocampus.numtrials(qdata)
    if end_trial === nothing
        end_trial = nt
    end
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    m_floor = (Hippocampus.floor_topology3(;nrefinements=3))
    poster_pos = Dict{Symbol, NTuple{3,Float64}}()
    # rearrnage so that the z-coordinate comes last
    for (k,v) in udata.poster_pos
        poster_pos[k] = v[[1,3,2]]
    end
    posters = Hippocampus.Posters(mm, poster_pos) 
    trialnr = Observable(start_trial)
    # poster name are sorted alpahbetically
    poster_ids = sort(collect(keys(poster_pos)))
    posteridx = udata.triggers[trialnr[],1] - 10
    poster_img = Observable(rotr90(load(Hippocampus.poster_img[poster_ids[posteridx]])))
    tg,gaze,pos, fixmask,fo = Hippocampus.get_trial(qdata,trialnr[];trial_start=1);
    gaze_point = Observable(Point3f(gaze[:,1]))
    pos_point = Observable(Point3f(pos[:,1]))
    line_of_gaze = Observable([Point3f(pos[:,1]), Point3f(gaze[:,1])])
    pos_spikes = Observable(Point3f[])
    gaze_spikes = Observable(Point3f[])
    i = Observable(1)
    on(trialnr) do _trialnr
        tg,gaze,pos, fixmask,fo = Hippocampus.get_trial(qdata,_trialnr;trial_start=1);
        posteridx = udata.triggers[_trialnr,1] - 10
        poster_img[] = rotr90(load(Hippocampus.poster_img[poster_ids[posteridx]]))
        i[] = 1
    end


    with_theme(Hippocampus.plot_theme) do
        fig = Figure()
        lscene = LScene(fig[1,1],show_axis=false)
        # creaete an inset axis with the current poster
        ax = Axis(fig[1,1], width=Relative(0.1), height=Relative(0.2), halign=0.1, valign=0.9, aspect=DataAspect())
        hidedecorations!(ax)
        image!(ax, poster_img)
        Hippocampus.plotmesh!(lscene, mm, color=:lightgray, showsegments=true, hide_floor=true, hide_ceiling=true)
        plot!(lscene, posters, shading=false)
        viz!(lscene, m_floor, color=:lightblue, showsegments=true)
        scatter!(lscene, pos_point)
        scatter!(lscene, gaze_point)
        scatter!(lscene, pos_spikes, color=:black)
        scatter!(lscene, gaze_spikes, color=:red)
        lines!(lscene, line_of_gaze)
        do_quit = Observable(false) 
        on(events(fig.scene).keyboardbutton) do event
            if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key == Keyboard.escape
                    do_quit[] = true
                end
            end
        end

        display(fig)
        record(fig.scene, fname;framerate=10) do io
            while true
                gaze_point[] = Point3f(gaze[:,i[]])
                pos_point[] = Point3f(pos[:,i[]])
                line_of_gaze[] = [Point3f(pos[:,i[]]), Point3f(gaze[:,i[]])]
                vidx = findfirst(vpvrp.placeviewidx[trialnr[]].==i[])
                if vidx !== nothing
                    if tg[i[]] > udata.timestamps[trialnr[],2] - udata.timestamps[trialnr[],1] # skip spikes during cue 
                        if fo[i[]] != "HintImage" # skip spikes during hint image
                            _pos_spikes = pos_spikes[]
                            push!(_pos_spikes, Point3f(pos[:,i[]]))
                            pos_spikes[] = _pos_spikes
                            _gaze_spikes = gaze_spikes[]
                            push!(_gaze_spikes, Point3f(gaze[:,i[]]))
                            gaze_spikes[] = _gaze_spikes
                        end
                    end
                end
                recordframe!(io)
                sleep(0.01)
                yield()
                i[] += 1
                if i[] > length(tg)
                    if trialnr[] < end_trial 
                        trialnr[] += 1
                    else
                        break
                    end
                    i[] = 1
                end
                if do_quit[]
                    break
                end
            end
        end
        fig
    end

end