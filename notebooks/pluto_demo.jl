### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 7e46e0b6-04c6-11eb-18f7-bf3077b15a72
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 2b54d586-0474-11eb-22bb-73ff01024985
begin
	using Revise
	using InteractiveUtils
	using PlutoUI
	using RvSpectMLBase
	using RvSpectML
	using EchelleInstruments
	using DataFrames
	using Query
	using Statistics
	using Plots
	using RvSpectMLPlots
end

# ╔═╡ a4d84168-048d-11eb-1c8f-0fa3136a1e0f
size(order_ccfs1,2)

# ╔═╡ f06adcc4-0480-11eb-08ba-5dc0b9d3062b
md"Mask width factor: $(@bind mask_scale_factor Slider(0.1:0.1:5; default=2, show_value=true))"

# ╔═╡ 8bf3f9bc-0479-11eb-144e-39bcfb215074
md"CCF Residuals? $(@bind plt_ccf_resid CheckBox(default=false))   
CCF Timeseries? $(@bind plt_ccf_timeseries CheckBox(default=false))   
CCF Orders? $(@bind plt_ccf_orders CheckBox(default=false))   
CCF Order Timeseries? $(@bind plt_ccf_order_timeseries CheckBox(default=false))"   

# ╔═╡ e0f9e742-048c-11eb-171f-4f1f3dea0f21
  begin
	if plt_ccf_orders
	local plt = plot()
		
	for ord in 1:size(order_ccfs1,2)
		if !plt_ccf_resid
			plot!(plt,v_grid_order_ccfs1, order_ccfs1[:,ord,1]./maximum(order_ccfs1[:,ord,1],dims=1) .+ 0.00*ord, label=:none)
			plot!(plt,v_grid_order_ccfs2, order_ccfs2[:,ord,1]./maximum(order_ccfs2[:,ord,1],dims=1) .+ 0.00*ord, label=:none)
		else
				plot!(plt,v_grid_order_ccfs1, order_ccfs1[:,ord,1]./maximum(order_ccfs1[:,ord,1],dims=1) .-  order_ccfs2[:,ord,1]./maximum(order_ccfs2[:,ord,1],dims=1) .+ 0.00*ord, label=:none)
		end
	end
    #local colorscale = cgrad(:balance)
    #plt = heatmap(v_grid_order_ccfs,collect(1:size(order_ccfs,2)),zvals', c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
    #add_time_gap_lines(plt,order_list_timeseries.times)
    xlabel!("v (m/s)")
    ylabel!("Order #")
    title!("CCF(v,t)-<CCF>(v) vs time")
    plt
	end
end

# ╔═╡ 4156a452-0476-11eb-2a93-6f6ba585ab3f
md"### Setup environment"

# ╔═╡ 37855416-0477-11eb-315c-4b4bca632ae5
begin
	expres_data_path = joinpath(homedir(),"Data/EXPRES/inputs")
	target_subdir = "101501"
 	fits_target_str = "101501"
 	output_dir = "examples/output/"
    df_files = EchelleInstruments.make_manifest(expres_data_path, target_subdir, EXPRES );
end;

# ╔═╡ 3ca4f348-0477-11eb-35a9-4ba6addfd1b4
begin
	# Configuration from param.in
	linelist_for_ccf_filename = "G8.espresso.mas"
	ccf_mid_velocity = -5e3
	max_spectra_to_use = 50
	df_files_use = df_files |>
	    @filter( _.target == fits_target_str ) |>
	    @take(max_spectra_to_use) |>
	    DataFrame;
end;

# ╔═╡ 8b7b2098-0486-11eb-2abe-c530be84ba4e
"MJD = $(df_files_use[obs_highlight,:bjd])  airmass = $(df_files_use[obs_highlight,:airmass])   epoch = $( df_files_use[obs_highlight,:expres_epoch])"

# ╔═╡ 4ad9e038-0477-11eb-2363-0f093d35f3ca
md"#### Read data files & extract each order"

# ╔═╡ 334533ae-0478-11eb-2c9e-07eedf379844
begin 
	pipeline_plan = PipelinePlan()
	dont_make_plot!(pipeline_plan, :movie)
end;

# ╔═╡ 4d93dcfc-0477-11eb-2619-c10f7fd1b909
begin
	spectra = map(row->EXPRES.read_data(row), eachrow(df_files_use))
	dont_need_to!(pipeline_plan,:read_spectra)
end;

# ╔═╡ abdf1064-0487-11eb-1f99-eb6563e505d7
md"Select Plot Obs# $(@bind obs_idx1 Slider(1:length(spectra),default=1,show_value=true))    Select Ref Obs# $(@bind obs_idx2 Slider(1:length(spectra),default=2,show_value=true))"

# ╔═╡ 5e945f0e-0490-11eb-1745-65666efed48b
obs_idx2

# ╔═╡ eb708bb0-0476-11eb-0ff3-bd45dc247e5a
function extract_orders(data::AbstractVector{ST}; verbose::Bool = false, orders_to_use = RvSpectMLBase.orders_to_use_default(first(data).inst) ) where { ST<:AbstractSpectra }
      order_list_timeseries = make_order_list_timeseries(data, orders_to_use=orders_to_use)
      order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
end;

# ╔═╡ 29638118-0475-11eb-19db-bdcac9172074
begin 
	olt1 = make_order_list_timeseries(spectra, orders_to_use=[obs_idx1])
    #olt1 = filter_bad_chunks(olt1,verbose=true)
	dont_need_to!(pipeline_plan,:extract_orders)
end

# ╔═╡ ff7bab4c-04be-11eb-3de1-27f14eef3b91
begin 
	olt2 = make_order_list_timeseries(spectra, orders_to_use=[obs_idx2])
    #olt1 = filter_bad_chunks(olt1,verbose=true)
end

# ╔═╡ 07c981e6-04bf-11eb-3c02-b9d8d6018719
begin 
	order_list_timeseries = make_order_list_timeseries(spectra, orders_to_use=[obs_idx1])
    #order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
end

# ╔═╡ b8eba050-04bf-11eb-16b7-358930a15ffd
RvSpectMLBase.make_orders_into_chunks(spectra[1], first(spectra).inst, orders_to_use = [5])


# ╔═╡ 059cc5ca-0478-11eb-3f29-6b50cb83c78e
line_list_df = prepare_line_list(linelist_for_ccf_filename, spectra, pipeline_plan,  v_center_to_avoid_tellurics=ccf_mid_velocity, Δv_to_avoid_tellurics = 30e3, recalc=true);

# ╔═╡ 08c3f4d2-047d-11eb-0b79-e1f9035e3df7
md"#### Actual calculations"

# ╔═╡ 0bae871e-0484-11eb-0b0e-234b276e02f3
begin
	lambda_min = minimum(line_list_df.lambda)
	lambda_max = maximum(line_list_df.lambda)
end;

# ╔═╡ 3c588554-0482-11eb-2fff-cd4645f903f5
@bind mask_lambda_min  Slider(lambda_min:10:lambda_max; default=lambda_min, show_value=true)

# ╔═╡ d93361c2-0483-11eb-3da2-db2a8ff59b90
@bind mask_lambda_max Slider(lambda_min:10:lambda_max; default=lambda_max, show_value=true) 

# ╔═╡ 17d6a8b4-0482-11eb-330b-cd9ed27943bf
line_list_use_df = @from i in line_list_df begin
    @where (mask_lambda_min <= i.lambda <= mask_lambda_max) 
	#@where (mask_weight_min <= i.weight <= mask_weight_max) 
    @select i
    @collect DataFrame
end;

# ╔═╡ 30571030-04c7-11eb-1516-dbd5589c12b9
function ccf_orders(order_list_timeseries::AbstractChunkListTimeseries,  line_list_df::DataFrame; verbose::Bool = false, range_no_mask_change::Real=30e3, ccf_mid_velocity::Real=0.0, mask_scale_factor::Real=1, mask_type::Symbol = :tophat )
      if mask_type == :tophat
         mask_shape = EchelleCCFs.TopHatCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      elseif mask_type == :gaussian
         mask_shape = EchelleCCFs.GaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :supergaussian
         mask_shape = EchelleCCFs.SuperGaussianCCFMask(order_list_timeseries.inst, σ_scale_factor=mask_scale_factor)
      elseif mask_type == :halfcos
         mask_shape = EchelleCCFs.CosCCFMask(order_list_timeseries.inst, scale_factor=mask_scale_factor)
      else
         @error("Requested mask shape (" * string(mask_type) * " not avaliable.")
      end

      line_list = EchelleCCFs.BasicLineList(line_list_df.lambda, line_list_df.weight)
      ccf_plan = EchelleCCFs.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list, midpoint=ccf_mid_velocity, range_no_mask_change=range_no_mask_change)
      v_grid = EchelleCCFs.calc_ccf_v_grid(ccf_plan)
      order_ccfs = EchelleCCFs.calc_order_ccf_chunklist_timeseries(order_list_timeseries, ccf_plan)
	return order_ccfs
end

# ╔═╡ f1def888-047b-11eb-1b72-63189bf2650a
(order_ccfs, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_df);

# ╔═╡ e22f9f90-04c7-11eb-2519-a74e022bcd38
order_ccfs

# ╔═╡ 85b04c2e-047c-11eb-330a-e321c5501b45
  begin
	if plt_ccf_orders
	local plt = plot()
	#local zvals = order_ccfs[:,:,obs_highlight]./maximum(order_ccfs[:,:,obs_highlight],dims=1).-0*plt_ccf_resid*mean(order_ccfs[:,:,obs_highlight]./maximum(order_ccfs[:,:,obs_highlight],dims=1),dims=2)
	local zvals = order_ccfs[:,:,obs_idx1]./maximum(order_ccfs[:,:,idx1],dims=1).- plt_ccf_resid.*reshape(mean(order_ccfs[:,:,obs_idx2]./maximum(order_ccfs[:,:,obs_idx2],dims=1),dims=3),size(order_ccfs)[1:2])
    local colorscale = cgrad(:balance)
    plt = heatmap(v_grid_order_ccfs,collect(1:size(order_ccfs,2)),zvals', c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
    add_time_gap_lines(plt,order_list_timeseries.times)
    xlabel!("v (m/s)")
    ylabel!("Order #")
    title!("CCF(v,t)-<CCF>(v) vs time")
    plt
	end
end

# ╔═╡ 5496cbaa-0476-11eb-0637-4d3d72087b10
md"### Pluto utils"

# ╔═╡ dd0b46c2-0472-11eb-08ef-b9b5984cf08c
@bind pluto_width Radio(["800px"=>"small", "1200px"=>"medium", "2400px"=>"large", "max-width"=>"max"], default="max-width")

# ╔═╡ d2a8c018-047e-11eb-2f3f-83b241d73c3d
begin
	mask_shape_keys = ["Gaussian", "Tophat", "Half Cos"]
	mask_shape_vals = [:gaussian, :tophat,:halfcos]
	mask_shape_dict = Dict(zip(mask_shape_keys,mask_shape_vals))
end;

# ╔═╡ 82a0f3b4-047f-11eb-2fb0-45fe12555aed
md"Mask shape: $(@bind mask_shape Radio(map(x->Pair(x[1],x[2]),zip(mask_shape_keys,mask_shape_vals)),default=\"Tophat\") )"

# ╔═╡ c623d82e-048c-11eb-3bde-43cb4234f60f
begin 
	(ccfs1, v_grid1) = ccf_total(olt1, line_list_use_df, pipeline_plan,  mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=false)
	#(order_ccfs1, v_grid_order_ccfs1) = ccf_orders(olt1, line_list_df,mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity)
end;

# ╔═╡ eeea64e8-0488-11eb-029e-8d0c711f4715
begin 
	local plt = plot() 
	plot!(plt, v_grid1,ccfs1./maximum(ccfs1,dims=1) .- 0.0 .* mean(ccfs1./maximum(ccfs1,dims=1),dims=2),labels=:none)
	#plot!(plt, v_grid_order_ccfs1, ccfs1./maximum(ccfs1,dims=1).-plt_ccf_resid.*ccfs2./maximum(ccfs2,dims=1) ,labels=:none)
	#plot!(plt, v_grid_order_ccfs2, ccfs2./maximum(ccfs2,dims=1).-plt_ccf_resid.*ccfs2./maximum(ccfs2,dims=1) ,labels=:none)
	#scatter!(plt, v_grid2, ccfs1[:,obs_idx1]./maximum(ccfs1[:,obs_idx1]) .- plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2),markersize=1.2,labels=:none)
	xlabel!("RV (m/s)")
	plt_ccf_resid ? ylabel!("CCF") : ylabel!("CCF")
	plt
end

# ╔═╡ 75c284dc-0489-11eb-34ff-d90b6385687a
ccfs1./maximum(ccfs1,dims=1)

# ╔═╡ 44c3e344-048e-11eb-2bac-cd263f61abae
begin 
	(ccfs2, v_grid2) = ccf_total(olt2, line_list_use_df, pipeline_plan,  mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=false)
	#(order_ccfs2, v_grid_order_ccfs2) = ccf_orders(olt2, line_list_df, pipeline_plan,mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity,recalc=true)
	#ccfs2 = reshape(sum(order_ccfs2,dims=3),(size(order_ccfs2,1),size(order_ccfs2,2)))
	#v_grid2 = v_grid_order_ccfs2
	#ccfs2
end;

# ╔═╡ ec169d3a-048e-11eb-178a-75dc8b1b3d95
 extrema(ccfs1.-ccfs2)

# ╔═╡ 5328c0d2-0478-11eb-357c-e750df71d6ff
begin 
	ENV["JULIA_NUM_THREADS"] = 4
	(ccfs, v_grid) = ccf_total(order_list_timeseries, line_list_use_df, pipeline_plan,  mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=false)
	#((ccfs, ccf_vars), v_grid) = ccf_total(order_list_timeseries, line_list_df, pipeline_plan,  mask_type=mask_shape_dict[mask_shape], mask_scale_factor=mask_scale_factor, ccf_mid_velocity=ccf_mid_velocity, recalc=true,  calc_ccf_var=true)
end;

# ╔═╡ 7f3e45ca-0478-11eb-2c41-4db39b42331e
begin 
	if plt_ccf_timeseries
		local plt = plot() 
	#plot!(plt, v_grid,ccfs./maximum(ccfs,dims=1) .- mean(ccfs./maximum(ccfs,dims=1),2),labels=:none)
	plot!(plt, v_grid, ccfs./maximum(ccfs,dims=1) .- plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2),labels=:none)
	scatter!(plt, v_grid, ccfs[:,obs_highlight]./maximum(ccfs[:,obs_highlight]) .- plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2),markersize=1.2,labels=:none)
	xlabel!("RV (m/s)")
	plt_ccf_resid ? ylabel!("CCF") : ylabel!("CCF")
	plt
	end
end

# ╔═╡ 78f39af8-0487-11eb-0948-bf2a888bd088
begin 
	if plt_ccf_orders
	local plt = plot() 
	#plot!(plt, v_grid,ccfs./maximum(ccfs,dims=1) .- mean(ccfs./maximum(ccfs,dims=1),2),labels=:none)
	plot!(plt, v_grid, ccfs./maximum(ccfs,dims=1) .- plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2),labels=:none)
	scatter!(plt, v_grid, ccfs[:,obs_highlight]./maximum(ccfs[:,obs_highlight]) .- plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2),markersize=1.2,labels=:none)
	xlabel!("RV (m/s)")
	plt_ccf_resid ? ylabel!("CCF") : ylabel!("CCF")
	plt
	end
end

# ╔═╡ 237fc58e-047c-11eb-0853-05ab5af54793
  begin
	if plt_ccf_timeseries
	local plt = plot()
	local zvals = ccfs./maximum(ccfs,dims=1).-plt_ccf_resid*mean(ccfs./maximum(ccfs,dims=1),dims=2)
    local colorscale = cgrad(:balance)
    plt = heatmap(v_grid,collect(1:size(ccfs,2)),zvals', c=colorscale, clims=(-maximum(abs.(zvals)),maximum(abs.(zvals))) )
    add_time_gap_lines(plt,order_list_timeseries.times)
    xlabel!("v (m/s)")
    ylabel!("Observation #")
    title!("CCF(v,t)-<CCF>(v) vs time")
    plt
	end
end

# ╔═╡ 2a705298-047a-11eb-3967-b3207243a5c3
begin
	mean_norm_ccf_total = mean(ccfs./maximum(ccfs,dims=1),dims=2)
	ccf1_norm = maximum(ccfs[:,obs_idx1])
end;

# ╔═╡ b2c72e50-0484-11eb-1d15-b953a0640a69
@bind mask_weight_min  Slider(0.0:0.05:0.95; default=0.0, show_value=true)

# ╔═╡ b27dcd6e-0484-11eb-0ab6-e73b9b4e80ec
@bind mask_weight_max Slider(0.05:0.05:1.0; default=1.0, show_value=true) 

# ╔═╡ a3737596-0473-11eb-0902-258a8ea5a515
HTML("<style> main { max-width:" * pluto_width * "; } </style> ")

# ╔═╡ 7ae40ec2-0475-11eb-3191-f75c848caaa0
HTML("""<style> input[type="range"] { width:250px; } </style>""")

# ╔═╡ 38a87a16-0475-11eb-1e64-b99500c57983
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ Cell order:
# ╠═c623d82e-048c-11eb-3bde-43cb4234f60f
# ╠═44c3e344-048e-11eb-2bac-cd263f61abae
# ╠═ec169d3a-048e-11eb-178a-75dc8b1b3d95
# ╠═5e945f0e-0490-11eb-1745-65666efed48b
# ╠═eeea64e8-0488-11eb-029e-8d0c711f4715
# ╠═e0f9e742-048c-11eb-171f-4f1f3dea0f21
# ╠═a4d84168-048d-11eb-1c8f-0fa3136a1e0f
# ╠═abdf1064-0487-11eb-1f99-eb6563e505d7
# ╠═f06adcc4-0480-11eb-08ba-5dc0b9d3062b
# ╠═82a0f3b4-047f-11eb-2fb0-45fe12555aed
# ╠═8bf3f9bc-0479-11eb-144e-39bcfb215074
# ╠═7f3e45ca-0478-11eb-2c41-4db39b42331e
# ╟─78f39af8-0487-11eb-0948-bf2a888bd088
# ╠═75c284dc-0489-11eb-34ff-d90b6385687a
# ╟─8b7b2098-0486-11eb-2abe-c530be84ba4e
# ╠═3c588554-0482-11eb-2fff-cd4645f903f5
# ╠═d93361c2-0483-11eb-3da2-db2a8ff59b90
# ╠═237fc58e-047c-11eb-0853-05ab5af54793
# ╟─4156a452-0476-11eb-2a93-6f6ba585ab3f
# ╠═7e46e0b6-04c6-11eb-18f7-bf3077b15a72
# ╠═2b54d586-0474-11eb-22bb-73ff01024985
# ╠═37855416-0477-11eb-315c-4b4bca632ae5
# ╠═3ca4f348-0477-11eb-35a9-4ba6addfd1b4
# ╟─4ad9e038-0477-11eb-2363-0f093d35f3ca
# ╠═334533ae-0478-11eb-2c9e-07eedf379844
# ╠═4d93dcfc-0477-11eb-2619-c10f7fd1b909
# ╠═eb708bb0-0476-11eb-0ff3-bd45dc247e5a
# ╠═29638118-0475-11eb-19db-bdcac9172074
# ╠═ff7bab4c-04be-11eb-3de1-27f14eef3b91
# ╠═07c981e6-04bf-11eb-3c02-b9d8d6018719
# ╠═b8eba050-04bf-11eb-16b7-358930a15ffd
# ╠═059cc5ca-0478-11eb-3f29-6b50cb83c78e
# ╟─17d6a8b4-0482-11eb-330b-cd9ed27943bf
# ╟─08c3f4d2-047d-11eb-0b79-e1f9035e3df7
# ╠═5328c0d2-0478-11eb-357c-e750df71d6ff
# ╠═0bae871e-0484-11eb-0b0e-234b276e02f3
# ╠═2a705298-047a-11eb-3967-b3207243a5c3
# ╠═30571030-04c7-11eb-1516-dbd5589c12b9
# ╠═e22f9f90-04c7-11eb-2519-a74e022bcd38
# ╠═f1def888-047b-11eb-1b72-63189bf2650a
# ╠═85b04c2e-047c-11eb-330a-e321c5501b45
# ╟─5496cbaa-0476-11eb-0637-4d3d72087b10
# ╠═dd0b46c2-0472-11eb-08ef-b9b5984cf08c
# ╠═d2a8c018-047e-11eb-2f3f-83b241d73c3d
# ╟─b2c72e50-0484-11eb-1d15-b953a0640a69
# ╟─b27dcd6e-0484-11eb-0ab6-e73b9b4e80ec
# ╟─a3737596-0473-11eb-0902-258a8ea5a515
# ╠═7ae40ec2-0475-11eb-3191-f75c848caaa0
# ╟─38a87a16-0475-11eb-1e64-b99500c57983
