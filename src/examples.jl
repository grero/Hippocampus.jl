using Hippocampus
using CairoMakie
using Meshes
using Colors
using LinearAlgebra


function theme_github_dark()
    _theme = theme_dark()
    _theme.backgroundcolor = "#0D1117"
    _theme.Axis.xgridvisible=false
    _theme.Axis.ygridvisible=false
    _theme
end

function theme_github_light()
    _theme = theme_light()
    _theme.Axis.xgridvisible=false
    _theme.Axis.ygridvisible=false
    _theme
end

function edge_smooth_example()
    m_floor = Hippocampus.repair_mesh(Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=2)));    
    # put a place field near one pillar
    occupancy = zeros(nelements(m_floor))
    weight = zeros(nelements(m_floor))

    point = (5.0, 8.0)
    # get the neighborhood
    kn = KNearestSearch(m_floor,7)
    idx = search(Meshes.Point(point...), kn) 

    occupancy = ones(nelements(m_floor))
    weight[idx] = rand(1:3, 7) 
    spm = Hippocampus.SpatialMapNew(m_floor, weight, occupancy)
    niter = 100
    spms = Hippocampus.SmoothedMap(spm;method=:laplace, α=0.1, niter=niter)
    # compare to a flat grid
    rg = RegularGrid((-12.5, -12.5), (12.5, 12.5);dims=(20,20))
    kn = KNearestSearch(rg,7)
    # find the same points in the grid as we used above
    idx_grid = [argmin([norm(p-pm) for pm in centroid.(rg)]) for p in centroid.(m_floor[idx])]
    occupancy_grid = ones(nelements(rg))
    weight_grid = zeros(nelements(rg))
    weight_grid[idx_grid] = weight[idx]
    spm_grid = Hippocampus.SpatialMapNew(rg, weight_grid, occupancy_grid) 
    spm_grids = Hippocampus.SmoothedMap(spm_grid;method=:laplace, α=0.1, niter=niter)
    bb = Hippocampus.find_boundary(m_floor)
    # gaussian smooth
    h = Meshes.ustrip(mean(sqrt.(measure.(mm))))
    D = [Meshes.ustrip(norm(p1-p2)) for p1 in centroid.(rg), p2 in centroid.(rg)]
    σ = h*sqrt((0.1*niter/2))
    @show σ
    F = exp.(-D.^2/(2*σ^2))
    Xs = F*weight_grid
    Ys = F*occupancy_grid
    Zs = Xs./Ys
    # 
    # overlay square mask
    r1 = Rect2f(2.5, 2.5, 5.0, 5.0)
    λs = Hippocampus.get_rate_map(spms)
    pidx = findall(λs .> eps())
    @show length(pidx)
    # find the correspondnce between each of the elements in the occluded floor to the unoccluded
    idx2 = [argmin([norm(p-pm) for pm in centroid.(rg)]) for p in centroid.(m_floor)]
    Z2 = zeros(nelements(rg))
    Z2[idx2] .= λs
    ΔZlg =  (Zs[idx2]-λs)
    ΔZlg2 = (Hippocampus.get_rate_map(spm_grids) - Zs)
    zmin,zmax = extrema(ΔZlg)
    midpoint = abs(zmin)/(zmax-zmin)
    zmin,zmax = extrema(ΔZlg2)
    midpoint2 = abs(zmin)/(zmax-zmin)

    segmentcolor = :gray75
    showsegments=false
    with_theme(theme_github_dark()) do
        fig = Figure(size=(656,309))
        ax1 = Axis(fig[1,1],aspect=1)
        Label(fig[1,1, TopLeft()], "A")
        hidedecorations!(ax1)
        hidespines!(ax1)
        v = Hippocampus.plotmap!(ax1, spm;showsegments=showsegments, segmentcolor=segmentcolor, colormap=:rain)
        Colorbar(fig[1,2], colorrange=extrema(v.color[]), colormap=:rain, label="Zo")
        ax2 = Axis(fig[1,3],aspect=1)
        Label(fig[1,3, TopLeft()], "B")
        hidedecorations!(ax2)
        hidespines!(ax2)
        Hippocampus.plotmap!(ax2, spms;showsegments=showsegments, segmentcolor=segmentcolor, colormap=:rain)
        Colorbar(fig[1,4], colorrange=extrema(λs), colormap=:rain, label="Zlo")
        ax3 = Axis(fig[2,1], aspect=1)
        Label(fig[2,1,TopLeft()], "D")
        hidedecorations!(ax3)
        hidespines!(ax3)
        Hippocampus.plotmap!(ax3, spm_grids;showsegments=showsegments, segmentcolor=segmentcolor, colormap=:rain)
        viz!(ax3, bb;color=:black,segmentsize=2.0)
        Colorbar(fig[2,2], colorrange=extrema(Hippocampus.get_rate_map(spm_grids)), colormap=:rain, label="Zlu")
        ax4 = Axis(fig[2,3], aspect=1)
        Label(fig[2,3, TopLeft()], "E")
        hidedecorations!(ax4)
        hidespines!(ax4)
        viz!(ax4, rg;color=Zs, colormap=:rain, showsegments=showsegments, segmentcolor=segmentcolor)
        viz!(ax4, bb;color=:black,segmentsize=2.0)
        Colorbar(fig[2,4], colorrange=extrema(Zs), colormap=:rain, label="Zgu")
        #poly!(ax4, r1, color=:black)
        #Hippocampus.plotmap!(ax4, spm_grids;showsegments=true, segmentcolor=:gray40, colormap=:rain)
        #show the difference map
        # TODO: We should show this on the occluded map, since that is where it matters in the end
        ax5 = Axis(fig[1,5], aspect=1)
        Label(fig[1,5,TopLeft()], "C")
        cm = colormap("RdBu", mid=midpoint)
        viz!(ax5, m_floor, color = ΔZlg, colormap=cm, showsegments=showsegments, segmentcolor=segmentcolor)
        viz!(ax5, bb;color=:black,segmentsize=2.0)
        Colorbar(fig[1,6], colormap=cm, label = "ΔZ")
        hidedecorations!(ax5)
        hidespines!(ax5)
        ax6 = Axis(fig[2,5], aspect=1)
        Label(fig[2,5, TopLeft()], "F")
        hidedecorations!(ax6)
        hidespines!(ax6)
        cm2 = colormap("RdBu", mid=midpoint2)
        viz!(ax6, rg, color=ΔZlg2, colormap=cm2)
        Colorbar(fig[2,6], colorrange=extrema(ΔZlg2), colormap=cm2, label="ΔZu")
        fig
    end
end