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

# ╔═╡ ba447bcc-0094-11eb-068b-6b5e045fd602
# Packages for basic plotting
begin 
	using Pkg
	Pkg.activate(".")
	#Pkg.develop("EchelleInstruments")
	#Pkg.instantiate()
	using RvSpectMLBase
	import EchelleInstruments
	import EchelleInstruments: EXPRES
	using DataFrames, Query
	using Statistics
	using Plots
	using PlutoUI
end

# ╔═╡ b29da99e-0142-11eb-3504-15cb11acd16f
# Packages for fitting GPs
begin 
	using LinearAlgebra
	using PDMats
	using StaticArrays
	using Stheno, TemporalGPs
	
	import Stheno: AbstractGP
	using Distributions
	import Distributions: AbstractMvNormal
	const Fx_PosteriorType = Distribution{Multivariate,Continuous}
end;

# ╔═╡ 78bc57ea-0140-11eb-2f53-15d0cb1e8b7b
md"# Visualize GP fit to spectra"

# ╔═╡ 19fc60d6-009a-11eb-3d9c-0f9152ece864
md"Divide by blaze: $(@bind div_blaze CheckBox(default=true))    Show tellurics? $(@bind plt_tellurics CheckBox())  Show continuum? $(@bind plt_continuum CheckBox())"
# Log y-axis? $(@bind log_y CheckBox()) # Need to fix GP fitting to logs

# ╔═╡ 26060010-0149-11eb-0863-33bfb7908253
md"Scaling factors for GP variance: $(@bind smooth_factor_var Slider(0.05:0.05:5; default=1.0, show_value=true))   & length scale $(@bind smooth_factor_length Slider(0.1:0.1:10; default=1.0, show_value=true))"

# ╔═╡ 9d4a5640-00a3-11eb-11f7-5d5e9c00fc61
md"### Boring code below"

# ╔═╡ ce23b456-00a5-11eb-0d07-2b15dfff6257
md"### Find files to read data from"

# ╔═╡ c5b8c394-013d-11eb-2353-7da0efded8f9
begin
	expres_data_path = joinpath(homedir(),"Data/EXPRES/inputs")
	target_subdir = "101501"
 	fits_target_str = "101501"
 	output_dir = "examples/output/"
    df_files = EchelleInstruments.make_manifest(expres_data_path, target_subdir, EXPRES );
end;

# ╔═╡ c364ee3e-0095-11eb-0faa-656ede64325f
begin
	# Configuration from param.in
	linelist_for_ccf_filename = "G8.espresso.mas"
	ccf_mid_velocity = -5e3
	max_spectra_to_use = 10
	df_files_use = df_files |>
	    @filter( _.target == fits_target_str ) |>
	    @take(max_spectra_to_use) |>
	    DataFrame;
end;

# ╔═╡ cfabe62c-0098-11eb-3b6d-ed7ab2f9f0f1
md"Select the observation index to plot: $(@bind obs_idx Slider(1:size(df_files_use,1); default=1, show_value=true))"

# ╔═╡ 99e1c70c-00a5-11eb-3cc9-fb8934931026
md"### Read data files & extract each order"

# ╔═╡ 01c2dfa6-0096-11eb-1361-af95574e6f95
begin 
	spectra = map(row->EXPRES.read_data(row), eachrow(df_files_use))
end;

# ╔═╡ 1a16d0ec-0095-11eb-3bac-0ddaaf387169
function extract_orders(data::AbstractVector{ST}; verbose::Bool = false, orders_to_use = RvSpectMLBase.orders_to_use_default(first(data).inst) ) where { ST<:AbstractSpectra }
      order_list_timeseries = make_order_list_timeseries(data, orders_to_use=orders_to_use)
      order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
end;

# ╔═╡ fdcf688a-0091-11eb-3fb4-9333baf67eab
 order_list_timeseries = extract_orders(spectra);

# ╔═╡ a6031a34-0098-11eb-38aa-891004881de9
md"Select a chunk ID: $(@bind chunk_idx Slider(1:num_chunks(order_list_timeseries); default=floor(Int,num_chunks(order_list_timeseries)//2), show_value=true))"

# ╔═╡ 4ab1c184-009b-11eb-3438-35502bf52add
begin 
	pixel_range_min = 1
	pixel_range_max = maximum(map(chid->length(order_list_timeseries[1][chid].λ),1:num_chunks(order_list_timeseries)))
	max_width_pixels = floor(Int,(pixel_range_max)//2)
	#pixel_range_max = min(length(order_list_timeseries[obs_idx][chunk_idx].λ), length(order_list_timeseries[obs_jdx][chunk_idx].λ))
	#max_width_pixels = floor(Int,minimum(map(chid->length(order_list_timeseries[1][chid].λ),1:num_chunks(order_list_timeseries)))//2)
	md"Central pixel: $(@bind mid_pixel Slider(pixel_range_min:pixel_range_max; default=floor(Int,(pixel_range_min+pixel_range_max)//2)))  and
	log2(Half-width in pixels): $(@bind log2_halfwidth_pixels Slider(0:0.2:log(max_width_pixels)/log(2); default=round(log(max_width_pixels)/log(2)*10)/10, show_value=true))"
end

# ╔═╡ 98c79090-014a-11eb-1ea9-c1b334406a3f
md"### Transform observed data"

# ╔═╡ 5b4656e8-014a-11eb-245f-65461ed5c02c
md"### Code for fitting GPs"

# ╔═╡ 5279c31c-0143-11eb-1e2c-bf94e52c0eb1
# Code to create prior GP
begin 
	#gp_param_default = [ 0.4890909216856761, 5.800274590507981e-5] .* smooth_factor   # TODO Generalize.  Values fr fit to one order of one EXPRES spectra
	#gp_param_default = [ 0.4890909216856761, 5.800274590507981e-5 * smooth_factor ]
	gp_param_default = [ 0.4890909216856761 * smooth_factor_var, 5.800274590507981e-5 * smooth_factor_length ]  
	gp_σ², gp_l = gp_param_default
	gp_kernel = gp_σ² * stretch(Matern52(), 1 / gp_l)
	prior_gp_naive = GP(gp_kernel, GPC())
	#f = to_sde(f_naive)   # if develop issues with StaticArrays could revert to this
	prior_gp = to_sde(prior_gp_naive, SArrayStorage(Float64))
end;

# ╔═╡ 6e3b4a28-0150-11eb-3f5c-59a21de3d298
# GP code assumes linear in flux.
log_y = false;

# ╔═╡ 44bce5e8-0147-11eb-18be-4df6dab237de
begin 
	# Make copy spectral chunks to plot, since many need to manipulate for plotting
	(pixels1, order1) = order_list_timeseries[obs_idx][chunk_idx].flux.indices
	flux1 = copy(order_list_timeseries[obs_idx][chunk_idx].flux)
	sigmasq_obs_trans = copy(order_list_timeseries[obs_idx][chunk_idx].var)
	# Either divide flux by blaze function or multiply continue by blaze
	blaze1 = order_list_timeseries.metadata[obs_idx][:blaze][pixels1,order1]
	if div_blaze
		flux1 ./= blaze1
		sigmasq_obs_trans ./= blaze1.^2
	end
	#flux1_norm = mean(flux1)
	#flux1 ./= flux1_norm
	#sigmasq_obs_trans ./= flux1_norm.^2 
	if log_y
		sigmasq_obs_trans ./= flux1.^2 
		flux1 = log10.(flux1)
	end
end;

# ╔═╡ d21e73e8-0147-11eb-0a9c-cd854964a3ee
# Compute GP posterior for entire order
begin
	flux1_norm = mean(flux1)
	flux_post = posterior(prior_gp(log.(order_list_timeseries[obs_idx][chunk_idx].λ), sigmasq_obs_trans./flux1_norm.^2), flux1./flux1_norm.-1)
end;

# ╔═╡ 4edfca32-0092-11eb-25d5-dbda1e798503
begin
	plt = plot(size=(800,500))
	if log_y  ylabel!("log Flux")   else ylabel!("Flux") end
	title!("Compare spectra to GP model")
	
	# Choose which pixels to plot
	halfwidth_pixels = floor(Int,2^log2_halfwidth_pixels)
	min_pixel = min(max(1,mid_pixel-halfwidth_pixels),length(order_list_timeseries[obs_idx][chunk_idx].flux))
	max_pixel = max(min(mid_pixel+halfwidth_pixels, length(order_list_timeseries[obs_idx][chunk_idx].flux)),1)

	xpred = order_list_timeseries[obs_idx][chunk_idx].λ
	pred_distr = marginals(flux_post(log.(xpred)))[min_pixel:max_pixel]
	mean_ypred = (mean.(pred_distr).+1.0)*flux1_norm
	sigma_ypred = sqrt.(var.(pred_distr)).*flux1_norm
	plot!(xpred[min_pixel:max_pixel], mean_ypred, color=:blue, label=:none )
	plot!(xpred[min_pixel:max_pixel], mean_ypred.+sigma_ypred, color=:cyan, linestyle=:dash, label=:none )
	plot!(xpred[min_pixel:max_pixel], mean_ypred.-sigma_ypred, color=:cyan, linestyle=:dash, label=:none )
	
	if plt_continuum
		continuum1 = copy(order_list_timeseries.metadata[obs_idx][:continuum][pixels1,order1])
		#continuum1 ./= flux1_norm
		if !div_blaze 
			continuum1 .*= blaze1
		end	
		if log_y
			continuum1 = log.(continuum1)
		end
		plot!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[min_pixel:max_pixel],continuum1[min_pixel:max_pixel], color=:grey, label=:none)
	end
	
	# Plot observations
	plot!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[min_pixel:max_pixel],flux1[min_pixel:max_pixel], color=:black, label=:none) # string(obs_idx))

	if plt_tellurics
		tellurics1 = order_list_timeseries.metadata[obs_idx][:tellurics][pixels1,order1]
		tellurics_idx1 = findall(x->x!=one(x),tellurics1[min_pixel:max_pixel]) .+ (min_pixel-1)
		scatter!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[tellurics_idx1],flux1[tellurics_idx1], markersize=1.2, color=:green, label=:none)
	end
	
	plt2 = plot(size=(800,300))
	plot!(xpred[min_pixel:max_pixel], flux1[min_pixel:max_pixel].-mean_ypred, color=:blue, label=:none )
	ylabel!("Residuals")
	xlabel!("λ (Å)")
	plot(plt,plt2,layout=(2,1), size=(800,800))
end

# ╔═╡ Cell order:
# ╟─78bc57ea-0140-11eb-2f53-15d0cb1e8b7b
# ╟─cfabe62c-0098-11eb-3b6d-ed7ab2f9f0f1
# ╟─a6031a34-0098-11eb-38aa-891004881de9
# ╟─4ab1c184-009b-11eb-3438-35502bf52add
# ╟─19fc60d6-009a-11eb-3d9c-0f9152ece864
# ╟─26060010-0149-11eb-0863-33bfb7908253
# ╠═4edfca32-0092-11eb-25d5-dbda1e798503
# ╟─9d4a5640-00a3-11eb-11f7-5d5e9c00fc61
# ╟─ba447bcc-0094-11eb-068b-6b5e045fd602
# ╟─ce23b456-00a5-11eb-0d07-2b15dfff6257
# ╠═c5b8c394-013d-11eb-2353-7da0efded8f9
# ╠═c364ee3e-0095-11eb-0faa-656ede64325f
# ╟─99e1c70c-00a5-11eb-3cc9-fb8934931026
# ╠═01c2dfa6-0096-11eb-1361-af95574e6f95
# ╠═1a16d0ec-0095-11eb-3bac-0ddaaf387169
# ╠═fdcf688a-0091-11eb-3fb4-9333baf67eab
# ╟─98c79090-014a-11eb-1ea9-c1b334406a3f
# ╠═44bce5e8-0147-11eb-18be-4df6dab237de
# ╟─5b4656e8-014a-11eb-245f-65461ed5c02c
# ╟─b29da99e-0142-11eb-3504-15cb11acd16f
# ╟─5279c31c-0143-11eb-1e2c-bf94e52c0eb1
# ╟─6e3b4a28-0150-11eb-3f5c-59a21de3d298
# ╟─d21e73e8-0147-11eb-0a9c-cd854964a3ee
