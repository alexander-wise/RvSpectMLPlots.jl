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
begin 
	using Pkg
	Pkg.activate(".")
	using RvSpectMLBase
	import EchelleInstruments
	import EchelleInstruments: EXPRES
	using DataFrames, Query
	using Statistics
	using Plots
	using PlutoUI
end

# ╔═╡ 78bc57ea-0140-11eb-2f53-15d0cb1e8b7b
md"# Interactive plots of 2 spectra"

# ╔═╡ 19fc60d6-009a-11eb-3d9c-0f9152ece864
md"""Divide by blaze: $(@bind div_blaze CheckBox())   Log y-axis? $(@bind log_y CheckBox())  Show tellurics? $(@bind plt_tellurics CheckBox())  Show continuum? $(@bind plt_continuum CheckBox())    $(@bind go_reset_range Button("Reset pixel range")) """

# ╔═╡ 53b7cab4-04e8-11eb-33c2-65008a5666a4
md"### Plotting functions"

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
	linelist_for_ccf_filename = "G8.espresso.mas"
	ccf_mid_velocity = -5e3
	max_spectra_to_use = 100
	df_files_use = df_files |>
	    @filter( _.target == fits_target_str ) |>
	    @take(max_spectra_to_use) |>
	    DataFrame;
end;

# ╔═╡ b00d8782-04f1-11eb-2c87-9d93fc5ec593
begin
	using FITSIO
	f = FITS(df_files_use.Filename[1])
	read_header(f[3])
end

# ╔═╡ 99e1c70c-00a5-11eb-3cc9-fb8934931026
md"### Read data files & extract each order"

# ╔═╡ 01c2dfa6-0096-11eb-1361-af95574e6f95
all_spectra = map(EXPRES.read_data,eachrow(df_files_use));

# ╔═╡ 1a16d0ec-0095-11eb-3bac-0ddaaf387169
function extract_orders(data::AbstractVector{ST}; verbose::Bool = false, orders_to_use = RvSpectMLBase.orders_to_use_default(first(data).inst) ) where { ST<:AbstractSpectra }
      order_list_timeseries = make_order_list_timeseries(data, orders_to_use=orders_to_use)
      order_list_timeseries = filter_bad_chunks(order_list_timeseries,verbose=true)
end;

# ╔═╡ fdcf688a-0091-11eb-3fb4-9333baf67eab
order_list_timeseries = extract_orders(all_spectra );

# ╔═╡ cfabe62c-0098-11eb-3b6d-ed7ab2f9f0f1
md"Select the observation indices to plot: $(@bind obs_idx Slider(1:length(order_list_timeseries); default=1, show_value=true))   $(@bind obs_jdx Slider(1:length(order_list_timeseries); default=2, show_value=true))"

# ╔═╡ a6031a34-0098-11eb-38aa-891004881de9
md"Select a chunk ID: $(@bind chunk_idx Slider(1:num_chunks(order_list_timeseries); default=floor(Int,num_chunks(order_list_timeseries)//2), show_value=true))   Order = $(order_list_timeseries[obs_idx].order[chunk_idx])"

# ╔═╡ 4ab1c184-009b-11eb-3438-35502bf52add
begin 
	go_reset_range
	pixel_range_min = 1
	pixel_range_max = maximum(map(chid->length(order_list_timeseries[1][chid].λ),1:num_chunks(order_list_timeseries)))
	max_width_pixels = floor(Int,(pixel_range_max)//2)
	md"Central pixel: $(@bind mid_pixel Slider(pixel_range_min:pixel_range_max; default=floor(Int,(pixel_range_min+pixel_range_max)//2)))  and
	Half-width in pixels: $(@bind halfwidth_pixels Slider(1:max_width_pixels; default=max_width_pixels))"
end

# ╔═╡ 066347f4-04d9-11eb-3a82-4f986133e7ac
get_meta(key::Symbol, obsid::Integer; digits::Integer = 3) = floor(parse(Float64,order_list_timeseries.metadata[obsid][key])*10^3)/10^3;

# ╔═╡ f226b23a-04d1-11eb-02ed-8d71509a5731
begin
	print_ti = round((order_list_timeseries.times[obs_idx]-order_list_timeseries.times[1])*100)/100
	print_tj = round((order_list_timeseries.times[obs_jdx]-order_list_timeseries.times[1])*100)/100
	pwv1 = get_meta(:PWV,obs_idx)
	pwv2 = get_meta(:PWV,obs_jdx)
	pwv2 = get_meta(:PWV,obs_jdx)
	air1 = get_meta(:AIRMASS,obs_idx)
	air2 = get_meta(:AIRMASS,obs_jdx)
	sun1 = get_meta(:SUNDIST,obs_idx)
	sun2 = get_meta(:SUNDIST,obs_jdx)
	moon1 = get_meta(:MOONDIST,obs_idx)
	moon2 = get_meta(:MOONDIST,obs_jdx)
	bc1 =  NaN # get_meta(:wtd_mdpt_bc,obs_idx)
	bc2 =  NaN # get_meta(:wtd_mdpt_bc,obs_jdx)
end;

# ╔═╡ f416a706-04e7-11eb-3fa1-53dc1bd0401d
string("Times = ", print_ti, " & ", print_tj, "   Δt = ", print_tj-print_ti, " days") 

# ╔═╡ 60ccca98-04f1-11eb-0b88-2993a6c13be7
string("BC  = ", bc1, " & ", bc2, "   ΔBC  = ", bc2-bc1)

# ╔═╡ f4d659c0-04e7-11eb-1f1a-05c3603a874d
string("PWV = ", pwv1, " & ", pwv2, "   ΔPWV = ", pwv2-pwv1) 

# ╔═╡ e0f528c8-04e7-11eb-135d-c19a28f13d58
string("Airmass = ", air1, " & ", air2, "   ΔAir = ", air2-air1) 

# ╔═╡ e60d0470-04e7-11eb-3457-7f4d35f8bffd
string("SunDist  = ", sun1, " & ", sun2, "   ΔD_Sun  = ", sun2-sun1, "    MoonDist = ", moon1, " & ", moon2, "   ΔD_Mood = ", moon2-moon1)

# ╔═╡ c1c3fb94-04cb-11eb-1830-f5a00c592959
md"### Some CSS tweaks to use more of screen"

# ╔═╡ eeb88092-015c-11eb-3989-2de8662a44bd

md"""Display width $(@bind pluto_width_pixels Select([ "600","800","1000","1200","1600","2000"] ))"""

# ╔═╡ 4edfca32-0092-11eb-25d5-dbda1e798503
function compare_two_spectra()
	plt = plot(size=(parse(Int,pluto_width_pixels),700))
	xlabel!("λ (Å)")
	ylabel!("Flux")
	title!("Comparing two spectra")
	# Make copy spectral chunks to plot, since many need to manipulate for plotting
	(pixels1, order1) = order_list_timeseries[obs_idx][chunk_idx].flux.indices
	(pixels2, order2) = order_list_timeseries[obs_jdx][chunk_idx].flux.indices
	flux1 = copy(order_list_timeseries[obs_idx][chunk_idx].flux)
	flux2 = copy(order_list_timeseries[obs_jdx][chunk_idx].flux)
	# Choose which pixels to plot
	min_pixel = min(max(1,mid_pixel-halfwidth_pixels),length(order_list_timeseries[obs_idx][chunk_idx].flux))
	max_pixel = max(min(mid_pixel+halfwidth_pixels, length(order_list_timeseries[obs_idx][chunk_idx].flux)),1)
	# Identify relevant blaze function and continuum
	blaze1 = order_list_timeseries.metadata[obs_idx][:blaze][pixels1,order1]
	blaze2 = order_list_timeseries.metadata[obs_jdx][:blaze][pixels2,order2]
	if plt_continuum
		continuum1 = copy(order_list_timeseries.metadata[obs_idx][:continuum][pixels1,order1])
		continuum2 = copy(order_list_timeseries.metadata[obs_jdx][:continuum][pixels2,order2])
	end
	# Either divide flux by blaze function or multiply continue by blaze
	if div_blaze
		flux1 ./= blaze1
		flux2 ./= blaze2
	else
		if plt_continuum
			continuum1 .*= blaze1
			continuum2 .*= blaze2
		end
	end
	# Calc scaling factor to make spectra line up
	scalefactor = mean(flux1)/mean(flux2)
	flux2 .*= scalefactor
	if plt_continuum
		continuum2 .*= scalefactor
	end
	if log_y
		flux1 = log10.(flux1)
		flux2 = log10.(flux2)
		if plt_continuum
			continuum1 = log10.(continuum1)
			continuum2 = log10.(continuum2)
		end
	end
	# Start plotting stuff
	plot!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[min_pixel:max_pixel],flux1[min_pixel:max_pixel], color=:blue, label=:none) # string(obs_idx))
	plot!(plt,order_list_timeseries[obs_jdx][chunk_idx].λ[min_pixel:max_pixel],flux2[min_pixel:max_pixel], color=:red, label=:none) # string(obs_jdx))
	if plt_tellurics
		tellurics1 = order_list_timeseries.metadata[obs_idx][:tellurics][pixels1,order1]
		tellurics_idx1 = findall(x->x!=one(x),tellurics1[min_pixel:max_pixel]) .+ (min_pixel-1)
		scatter!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[tellurics_idx1],flux1[tellurics_idx1], markersize=1.2, color=:green, label=:none)
		tellurics2 = order_list_timeseries.metadata[obs_idx][:tellurics][pixels2,order2]
		tellurics_idx2 = findall(x->x!=one(x),tellurics2[min_pixel:max_pixel]) .+ (min_pixel-1)
		scatter!(plt,order_list_timeseries[obs_jdx][chunk_idx].λ[tellurics_idx2],flux2[tellurics_idx2], markersize=1.2, color=:green, label=:none)
	end
	if plt_continuum
		plot!(plt,order_list_timeseries[obs_idx][chunk_idx].λ[min_pixel:max_pixel],continuum1[min_pixel:max_pixel], color=:cyan, label=:none)
		plot!(plt,order_list_timeseries[obs_jdx][chunk_idx].λ[min_pixel:max_pixel],continuum2[min_pixel:max_pixel], color=:magenta, label=:none)
	end
	#=
	if plt_mask
		for x in lines_in_mask
			plot([x,x],[ylims(plt)[1], ylims(plt)[2]],color=:black, label=:none)
		end
	end
	=#
	plt
end;

# ╔═╡ b6a6ad72-04d7-11eb-0bd5-81c1ccc85a5c
compare_two_spectra()

# ╔═╡ 9690e8b2-04d8-11eb-1449-a78cf2d8f32f
HTML("<style> main { max-width:" * string(pluto_width_pixels) * "px; } </style> ")

# ╔═╡ d22660c6-015b-11eb-0f6d-8745e4cd4c35
HTML("""<style> input[type="range"] { width:240px; } </style>""")

# ╔═╡ Cell order:
# ╟─78bc57ea-0140-11eb-2f53-15d0cb1e8b7b
# ╟─cfabe62c-0098-11eb-3b6d-ed7ab2f9f0f1
# ╟─a6031a34-0098-11eb-38aa-891004881de9
# ╟─4ab1c184-009b-11eb-3438-35502bf52add
# ╟─19fc60d6-009a-11eb-3d9c-0f9152ece864
# ╟─b6a6ad72-04d7-11eb-0bd5-81c1ccc85a5c
# ╟─f416a706-04e7-11eb-3fa1-53dc1bd0401d
# ╟─60ccca98-04f1-11eb-0b88-2993a6c13be7
# ╟─f4d659c0-04e7-11eb-1f1a-05c3603a874d
# ╟─e0f528c8-04e7-11eb-135d-c19a28f13d58
# ╟─e60d0470-04e7-11eb-3457-7f4d35f8bffd
# ╟─53b7cab4-04e8-11eb-33c2-65008a5666a4
# ╟─4edfca32-0092-11eb-25d5-dbda1e798503
# ╟─066347f4-04d9-11eb-3a82-4f986133e7ac
# ╟─f226b23a-04d1-11eb-02ed-8d71509a5731
# ╟─9d4a5640-00a3-11eb-11f7-5d5e9c00fc61
# ╟─ba447bcc-0094-11eb-068b-6b5e045fd602
# ╟─ce23b456-00a5-11eb-0d07-2b15dfff6257
# ╟─c5b8c394-013d-11eb-2353-7da0efded8f9
# ╟─c364ee3e-0095-11eb-0faa-656ede64325f
# ╟─99e1c70c-00a5-11eb-3cc9-fb8934931026
# ╟─01c2dfa6-0096-11eb-1361-af95574e6f95
# ╟─1a16d0ec-0095-11eb-3bac-0ddaaf387169
# ╠═fdcf688a-0091-11eb-3fb4-9333baf67eab
# ╠═b00d8782-04f1-11eb-2c87-9d93fc5ec593
# ╟─c1c3fb94-04cb-11eb-1830-f5a00c592959
# ╟─eeb88092-015c-11eb-3989-2de8662a44bd
# ╟─9690e8b2-04d8-11eb-1449-a78cf2d8f32f
# ╟─d22660c6-015b-11eb-0f6d-8745e4cd4c35
