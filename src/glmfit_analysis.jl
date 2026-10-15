using Hippocampus
using CairoMakie
using Meshes
using CairoMakie.Colors
using StatsBase

function plot_theme()
    theme = theme_minimal()
    theme.Axis.xticksvisible = true
    theme.Axis.yticksvisible = true
    theme.Colorbar.ticksvisible = true
    theme
end

function hist_diagonal!(ax, x::AbstractVector{<:Real}, y::AbstractVector{<:Real};offset=0)
    z = y-x
    hh = fit(Histogram, z)
    s = 0.9*step(hh.edges[1])
    # Draw rotated histogram bars manually
    θ = -π/4
    rot_matrix = [cos(θ) -sin(θ); sin(θ) cos(θ)] 
    v = Vec2(offset*[cos(-θ),sin(-θ)])
    @show offset
    for (bin_edge, height) in zip(hh.edges[1][1:end-1], 20*hh.weights)
        # Calculate rotated bar coordinates
        x1, y1 = rot_matrix*[bin_edge, offset]
        x2, y2 = rot_matrix*[bin_edge+s,offset]
        x3,y3 = rot_matrix*[bin_edge+s, offset+height]
        x4,y4 = rot_matrix*[bin_edge, offset+height]
        poly!(ax, [x1, x2, x3, x4], [y1, y2, y3, y4],
            color = :blue, strokecolor = :black)
    end
    # line
    x1,y1 = rot_matrix*[hh.edges[1][1] - 0.5*s, offset-2.5]
    x2,y2 = rot_matrix*[hh.edges[1][end] + 0.5*s, offset-2.5]
    lines!(ax, [x1,x2],[y1,y2], color=:black)
    # get the tick locations
    tl,xmin, xmax = Makie.optimize_ticks(extrema(z)...)
    # draw the ticks
    tl_xy0 = Point2f.(eachcol(rot_matrix*permutedims([tl fill(offset-2.5, length(tl))])))
    tl_xy1 = Point2f.(eachcol(rot_matrix*permutedims([tl fill(offset-10.0, length(tl))])))
    linesegments!(ax, collect(zip(tl_xy0, tl_xy1)),color=:black)
    #text!(ax, tl_xy1, text=string.(tl), align=(:center, :top),rotation=-π/4)
end

function hist_diagonal(x::AbstractVector{<:Real}, y::AbstractVector{<:Real};kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1],aspect=1)
    hist_diagonal!(ax, x, y;kwargs...)
    fig
end

function plot_rf_analysis!(lg, celldir;nrefinements=(p=2,g=2))
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=nrefinements.p))
    rf_spatial,jm = cd(celldir) do
        rf = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 1000;nrefinements=nrefinements,smooth=true, smoothing_method=:laplace, α=0.1, niter=50, redo=fname->false, min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02,trial_start=2, pv_threshold=0.001)
        jm = Hippocampus.JointMap(;redo=fname->false, do_save=true, nrefinements=nrefinements,min_speed=1.0, min_place_obs=5, min_view_obs=5, min_place_duration=0.05, min_view_duration=0.02, trial_start=2)
        rf,jm
    end
    spm = Hippocampus.SpatialMapNew(jm,m_floor)
    #plot the occupancy, spike counts, raw firing rates, and the smoothed firing rate    
    lg1 = GridLayout(lg[1,1])
    Hippocampus.plot_response_fields!(lg1, rf_spatial, spm.occupancy;colormap=:viridis, label="Occupancy [s]", colorbar_below=true)
    lg2 = GridLayout(lg[1,2])
    Hippocampus.plot_response_fields!(lg2, rf_spatial, spm.weight;colormap=:viridis, label="Spike count", colorbar_below=true)
    lg3 = GridLayout(lg[1,3])
    Hippocampus.plot_response_fields!(lg3, rf_spatial, Hippocampus.get_rate_map(spm);colormap=:viridis, label="Raw spike rate [Hz]", colorbar_below=true)
    lg4 = GridLayout(lg[1,4])
    Hippocampus.plot_response_fields!(lg4, rf_spatial;colormap=:viridis, label="Smooth spike rate [Hz]", colorbar_below=true)

end

function plot_glmfit_new!(lg, ll, ll0, β,nrefinements,α;hidx=1,non_visited::AbstractVector{Int64}=Int64[])
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=nrefinements.p))
    # out-data log-likelihood
    lg2 = GridLayout(lg[1,1])
    ax1 = Axis(lg2[1,1],xticks=WilkinsonTicks(3))
    color = fill(parse(Colorant,:black), length(ll0))
    color[hidx] = parse(Colorant, :red)
    scatter!(ax1, -ll0, -ll, color=color)
    ablines!(ax1, 0.0, 1.0, color=:black, linestyle=:dot)
    # add an axis showing a histogram
    #offset = -(maximum(maximum.([-ll0,-ll])) - minimum(minimum.([-ll0, -ll])))
    #offset = maximum(maximum.([-ll0,-ll]))
    #offset = -0.97*sqrt(maximum(-ll0)^2 + maximum(-ll)^2)
    #hist_diagonal!(ax1, -ll,-ll0;offset=offset)
    ax1.xlabel = "LL0"
    ax1.ylabel = "LL"
    ax3 = Axis(lg2[1,2], xticks=WilkinsonTicks(3))
    hist!(ax3, -ll+ll0, color=:gray)
    ax3.xlabel = "ll - ll0"

    lg1 = GridLayout(lg[1,2])
    ax2 = Axis(lg1[1,1], aspect=1)
    hidedecorations!(ax2)
    ax2.bottomspinevisible = false
    ax2.leftspinevisible = false
    ax2.xticklabelsvisible = false
    ax2.yticklabelsvisible = false
    viz!(ax2, m_floor;color=:lightgray)
    zz = β[2:end,hidx]
    zz[non_visited] .= NaN
    viz!(ax2, m_floor;color=zz)
    Colorbar(lg1[2,1], colorrange=extrema(filter(isfinite, zz)),vertical=false, flipaxis=false, ticks=WilkinsonTicks(3),label="β",ticksvisible=true)
    Label(lg[1,3], "α = $(α)", rotation=-π/2, tellwidth=true, tellheight=false)
    colsize!(lg, 2, Relative(0.25))
end

function plot_glmfit_new(celldir::String;figsize=(600,1000), nrefinements=(p=2,g=2), α=[0.1, 1.0, 10.0, 100.0, 1000.0],nruns=20, kwargs...)
    sessiondir = Hippocampus.DPHT.get_level_path("session", celldir)
    jocc = cd(sessiondir) do
        Hippocampus.JointOccupancy(;redo=fname->false, do_save=true, nrefinements=nrefinements,trial_start=2,min_speed=-1.0)
    end
    Z = Hippocampus.get_spatial_map(jocc, nrefinements)
    non_visited = findall(Z.==0.0)
    ll,ll0, β = Hippocampus.run_glm_zip(celldir, nrefinements,α;nruns=nruns)
    with_theme(plot_theme()) do
        fig = Figure(size=figsize)
        lg0 = GridLayout(fig[1,1])
        plot_rf_analysis!(lg0, celldir;nrefinements=nrefinements)
        lgs = [GridLayout(fig[i+1,1]) for i in 1:length(α)]
        for (ii,lg) in enumerate(lgs)
            plot_glmfit_new!(lg, ll[:,ii], ll0[:,ii], β[:,:,ii], nrefinements,α[ii];non_visited=non_visited,kwargs...)
        end
        fig
    end
end

function plot_glmfit_new(ll, ll0, β,nrefinements;figsize=(700,400))
    with_theme(theme_minimal()) do
        fig = Figure(size=figsize)
        lg = GridLayout(fig[1,1])
        plot_glmfit_new!(lg, ll, ll0, β, nrefinements)
        fig
    end
end

