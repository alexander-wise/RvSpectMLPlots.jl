module RvSpectMLPlots

using RvSpectMLBase
using Plots
using ColorSchemes

using LinearAlgebra
using StatsBase

include("spectra.jl")
export plot_spectrum_chunks, get_λs, add_time_gap_lines

include("scalpels.jl")
export make_plots_scalpels

include("dcpca.jl")
export plot_basis_vectors, plot_basis_scores, plot_basis_scores_cor

end
