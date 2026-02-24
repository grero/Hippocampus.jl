module PaperFigures
using CairoMakie
using Meshes
using Hippocampus

plot_theme = Theme(Axis=(xlabelsize=14, ylabelsize=14,
                           xticklabelsize=14, yticklabelsize=14,
                           topspinevisible=false, rightspinevisible=false,
                           xgridvisible=false, ygridvisible=false,ylabelvisible=true,
                           xticklabelsvisible=true, xlabelvisible=true),
                     Scatter=(markersize=10px,),
                     Lines=(linewidth=3,),
                     fontsize=14)

"""
Place cells
"""
function figure2(spatial_cells::Vector{String},example_idx::Vector{Int64})
    m_floor = Shadow("xy")(Hippocampus.floor_topology3(;nrefinements=3))
    with_theme(plot_theme) do
        fig = Figure(size=(900,800))
        # Example of a cell with place activity
        lgm = GridLayout(fig[1,1])
        for ii in example_idx
            lg1 = GridLayout(lgm[1,ii])
            spm,sic,rfs = cd(spatial_cells[ii]) do
                spm = Hippocampus.SpatialMapNew()
                sic = Hippocampus.compute_skaggs_sic(Hippocampus.SpatialInformationContent, 10_000;load_only=true, smooth=true, smoothing_method=:laplace, σ=0.1, niter=100)
                rfs = Hippocampus.get_response_fields(Hippocampus.SpatialResponseFields, 10_000;smooth=true, smoothing_method=:laplace, α=0.1, niter=100,pv_threshold=0.01) 
                spm,sic,rfs
            end
            sml = Hippocampus.SmoothedMap(spm;smooth_method=:laplace,α=0.1, niter=100)
            lg12 = GridLayout(lg1[1,1])
            ax = Axis(lg12[1,1],aspect=1)
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
            Colorbar(lg12[2,1], colorrange=extrema(filter(isfinite, Z)), label="Firing rate [Hz]",vertical=false,tellwidth=false)
            #colsize!(lg1, 1, Relative(0.6))

            # show SIC distribution
            #lg2 = GridLayout(fig[1,2])
            ax2 = Axis(lg12[3,1])
            Label(lg12[3,1, TopLeft()],"B")
            boxplot!(ax2, fill(1.0, length(sic.sic)), sic.sic,color=:gray,show_notch=true,orientation=:horizontal,show_outliers=false)
            vlines!(ax2, sic.sic0, linestyle=:dot, color=:black)
            ax2.yticklabelsvisible = false
            ax2.yticksvisible = false
            ax2.bottomspinevisible = true 
            ax2.leftspinevisible = false
            ax2.rightspinevisible = false 
            ax2.xlabel = "SIC"
            rowsize!(lg12, 1, Relative(0.8))
        end

        # summary showing total number of fields and coverage
        lg3 = GridLayout(fig[2,1]) 
        Hippocampus.plot_n_fields!(lg3, Hippocampus.SpatialResponseFields,spatial_cells;labels=["C","D","E"], pv_threshold=0.01, smooth=true, smoothing_method=:laplace, α=0.1,niter=100)
        rowsize!(fig.layout, 1, Relative(0.6))
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

end #module