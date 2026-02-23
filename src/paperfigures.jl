using CairoMakie
using Meshes


"""
Place cells
"""
function figure2(spatial_cells::Vector{String})
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3))
    with_theme(Hippocampus.plot_theme) do
        fig = Figure(size=(800,700))
        # Example of a cell with place activity
        lg1 = GridLayout(fig[1,1])
        spm,sic,rfs = cd(allcelldirs[39]) do
            spm = Hippocampus.SpatialMapNew()
            sic = Hippocampus.compute_skaggs_sic(Hippocampus.SpatialInformationContent, 10_000;load_only=true)
            rfs = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=100,pv_threshold=0.01) 
            spm,sic,rfs
        end
        sml = Hippocampus.SmoothedMap(spm;smooth_method=:laplace,α=0.1, niter=100)
        ax = Axis(lg1[1,1],aspect=1)
        Label(lg1[1,1,TopLeft()], "A")
        hidedecorations!(ax)
        ax.leftspinevisible = false
        ax.bottomspinevisible = false
        # show the outline of the maze
        viz!(ax, m_floor;color=:lightgray)
        Z = Hippocampus.get_rate_map(sml)
        viz!(ax,m_floor;color=Z)

        clusters = Hippocampus.merge_fields(rfs)
        nclusters = Hippocampus.get_num_fields(rfs)
        # relative number of times we get cluster of at least length.(clusters) randomly
        threshold = dropdims(sum(nclusters,dims=2),dims=2)/size(nclusters,2)
        # only keep fields where the probabilty of getting the same field in the surroages is less than 0.01
        valid_cluster_idx = findall(threshold .< 0.01)
        for cluster in clusters[valid_cluster_idx] 
            bb = Hippocampus.find_boundary(m_floor, rfs.binidx[cluster])
            viz!(ax, bb;color=:black)
        end
        Colorbar(lg1[1,2], colorrange=extrema(filter(isfinite, Z)), label="Firing rate [Hz]")

        # show SIC distribution
        lg2 = GridLayout(fig[1,2])
        ax2 = Axis(lg2[1,1],yaxisposition=:right)
        Label(lg2[1,1, TopLeft()],"B")
        boxplot!(ax2, fill(1.0, length(sic.sic)), sic.sic,color=:gray,show_notch=true)
        hlines!(ax2, sic.sic0, linestyle=:dot, color=:black)
        ax2.xticklabelsvisible = false
        ax2.xticksvisible = false
        ax2.bottomspinevisible = false
        ax2.ylabel = "SIC"
        colsize!(fig.layout, 2, Relative(0.2))

        # summary showing total number of fields and coverage
        lg3 = GridLayout(fig[2,1:2]) 
        Hippocampus.plot_n_fields!(lg3, Hippocampus.SpatialResponseFields,spatial_cells;labels=["C","D","E"], pv_threshold=0.01, smooth=true, smoothing_method=:laplace, α=0.1,niter=100)
        fig
    end
end


"""
View cells
"""
function figure3()
end


"""
Mixed selective and conjunction cells
"""
function figure4()
end