"""
Convenience code for plotting CCFs

Author: Eric Ford
Date:   October 2020
"""

"""   `make_plot_ccf_vs_time(ccfs; v_grid, t_idx, output_path, save_fig )`
Overplot CCFs (2d array) at multiple times.
Inputs:
- `ccfs`:  2d array of CCFs at different times
Optional Inputs:
- `v_grid`: velocities that CCF was evaluated at
- `t_idx`:  restrict plotting to these time indicies
- `save_fig`: save figures to file if true
- `output_path`: prefix to filename
Returns:
- plot object
"""
function make_plot_ccf_vs_time(ccfs::AbstractArray{T2,2}
                ; v_grid::AbstractRange = 1:size(ccfs,1),
                t_idx = 1:size(ccfs,2),
                output_path::String = "", save_fig::Bool = false ) where
                { T1<:Real, T2<:Real, T3<:Real, T4<:Real  }

    @assert length(v_grid) == size(ccfs,1)
    @assert length(t_idx) >= 1

    plt = plot()
    plot!(plt,v_grid,ccfs[:,t_idx]./maximum(ccfs[:,t_idx],dims=1),label=:none)
    xlabel!("v (m/s)")
    ylabel!("CCF")
    if save_fig   savefig(plt,output_path * "_ccf_sum.png")   end
    return plt
end


"""   `make_heatmap_ccf_vs_time(ccfs; v_grid, t_idx, output_path, save_fig )`
Plots heatmap of CCFs at multiple times minus the mean CCF (averaged over time).
Inputs:
- `ccfs`:  2d array of CCFs at different times
Optional Inputs:
- `v_grid`: velocities that CCF was evaluated at
- `save_fig`: save figures to file if true
- `output_path`: prefix to filename
Returns:
- plot object
"""
function make_heatmap_ccf_vs_time(ccfs::AbstractArray{T2,2}
                ; v_grid::AbstractRange = 1:size(ccfs,1),
                times::AbstractVector{T4} = collect(1:length(rvs)),
                output_path::String = "", save_fig::Bool = false ) where
                { T1<:Real, T2<:Real, T3<:Real, T4<:Real  }

    @assert length(v_grid) == size(ccfs,1)

    ccfs_minus_mean = ccfs./maximum(ccfs,dims=1) .- mean( ccfs./ maximum(ccfs,dims=1),dims=2)
    zvals = ccfs_minus_mean
    colorscale = cgrad(:balance)
    plt = plot()
    heatmap!(plt, v_grid, 1:size(ccfs,2),zvals',c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
    add_time_gap_lines(plt,times)
    title!("CCF(v,t)-<CCF> vs Time")
    xlabel!("v (m/s)")
    ylabel!("Observation #")
    if save_fig
        savefig(joinpath(output_path,"_ccf_sum_vs_time_heatmap.png"))
    end
    return plt
end



"""   `make_heatmap_ccf_vs_order(ccfs; v_grid, t_idx, output_path, save_fig )`
Plots heatmap of CCFs for multiple overs minus the mean CCF (averaged over time).
Inputs:
- `ccfs`:  2d array of CCFs at different times
Optional Inputs:
- `v_grid`: velocities that CCF was evaluated at
- `save_fig`: save figures to file if true
- `output_path`: prefix to filename
Returns:
- plot object
"""
function make_heatmap_ccf_vs_order(order_ccfs::AbstractArray{T2,3}
                ; v_grid::AbstractRange = 1:size(order_ccfs,1),
                #= times::AbstractVector{T4} = collect(1:length(rvs)), =#
                order_labels::AbstractVector{T5} = 1:size(order_ccfs,2),
                output_path::String = "", save_fig::Bool = false ) where
                { T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real  }

    @assert length(v_grid) == size(order_ccfs,1)

    obs = 1:size(order_ccfs,2)
    #order_labels = map(c->order_list_timeseries.chunk_list[1].data[c].λ.indices[2], 1:size(order_ccfs,2))
    orders_to_plot = findall(c->sum(order_ccfs[:,c,obs])>0, 1:size(order_ccfs,2))
    zvals =  reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=3)./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=3),dims=1),size(order_ccfs,1),size(order_ccfs[:,orders_to_plot,obs],2)) .-
             reshape(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))./maximum(sum(order_ccfs[:,orders_to_plot,obs],dims=(2,3))),size(order_ccfs,1) )
    colorscale = cgrad(:balance)
    plt = plot()
    heatmap!(plt,v_grid,order_labels[orders_to_plot], zvals', clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
    #println("size(v_grid) = ", size(v_grid_order_ccfs), "  size(order_labels[plt]) = ", size(order_labels[orders_to_plot]), "  size(zvals') = ", size(zvals') )
    xlabel!("v (m/s)")
    ylabel!("Order ID")
    title!("CCF-<CCF> for order_idx = " * string(order_labels[minimum(orders_to_plot)]) * " - " * string(order_labels[maximum(orders_to_plot)]))
    if save_fig
            if save_plot(pipeline_plan,:ccf_orders)   savefig(plt,output_path * " _ccf_orders.png")   end
    end

    return plt
end
