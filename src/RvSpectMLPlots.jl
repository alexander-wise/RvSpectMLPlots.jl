module RvSpectMLPlots

using RvSpectMLBase
using Plots
using ColorSchemes

using LinearAlgebra
using StatsBase

include("spectra.jl")
export plot_spectrum_chunks, get_λs, add_time_gap_lines

include("ccf.jl")
export make_plot_ccf_vs_time, make_heatmap_ccf_vs_time
#export make_plot_ccf_vs_order, make_plot_ccf_vs_chunk  # TODO: Write
export make_heatmap_ccf_vs_order

include("scalpels_code.jl")
export make_plots_scalpels

include("dcpca.jl")
export plot_basis_vectors, plot_basis_scores, plot_basis_scores_cor

end
