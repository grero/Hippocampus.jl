using Makie
using Meshes
using FileIO

function visualize_trial(qdata::Hippocampus.UnityRaytraceData, udata::Hippocampus.UnityData)
    mm = Hippocampus.get_maze_mesh(;nrefinements=2)
    m_floor = (Hippocampus.floor_topology3(;nrefinements=3))
    poster_pos = Dict{Symbol, NTuple{3,Float64}}()
    # rearrnage so that the z-coordinate comes last
    for (k,v) in udata.poster_pos
        poster_pos[k] = v[[1,3,2]]
    end
    posters = Hippocampus.Posters(mm, poster_pos) 
    # poster name are sorted alpahbetically
    poster_ids = sort(collect(keys(poster_pos)))
    posteridx = udata.triggers[1,1] - 10
    poster_img = load(Hippocampus.poster_img[poster_ids[posteridx]])
    tg,gaze,pos, fixmask,fo = Hippocampus.get_trial(qdata,1;trial_start=1);
    gaze_point = Observable(Point3f(gaze[:,1]))
    pos_point = Observable(Point3f(pos[:,1]))
    line_of_gaze = Observable([Point3f(pos[:,1]), Point3f(gaze[:,1])])

    with_theme(Hippocampus.plot_theme) do
        fig = Figure()
        lscene = LScene(fig[1,1],show_axis=false)
        # creaete an inset axis with the current poster
        ax = Axis(fig[1,1], width=Relative(0.1), height=Relative(0.2), halign=0.1, valign=0.9, aspect=DataAspect())
        hidedecorations!(ax)
        image!(ax, rotr90(poster_img))
        Hippocampus.plotmesh!(lscene, mm, color=:lightgray, showsegments=true, hide_floor=true, hide_ceiling=true)
        plot!(lscene, posters, shading=false)
        viz!(lscene, m_floor, color=:lightblue, showsegments=true)
        scatter!(lscene, pos_point)
        scatter!(lscene, gaze_point)
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
        i = 1
        while true
            gaze_point[] = Point3f(gaze[:,i])
            pos_point[] = Point3f(pos[:,i])
            line_of_gaze[] = [Point3f(pos[:,i]), Point3f(gaze[:,i])]
            sleep(0.01)
            yield()
            i += 1
            if i > length(tg)
                i = 1
            end
            if do_quit[]
                break
            end
        end
        fig
    end

end