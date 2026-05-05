module LLRParsing

using HiRepParsing
using Parsers
using HDF5
using StatsBase
using Statistics
using LinearAlgebra
using MadrasSokal
using ProgressMeter
using NaturalSort
using Parsers
using Plots
using LaTeXStrings
using Peaks
using Roots
using LsqFit
using PCHIPInterpolation
# compression support
using TranscodingStreams
using CodecZstd

include("errors.jl")
export mean_std_of_mean_cov, mean_std_of_mean, apply_jackknife
include("errorstring.jl")
export errorstring
include("parse_std.jl")
export importance_sampling_dir_hdf5
include("std_observables.jl")
export std_observables
include("parse_llr.jl")
export get_repeat_and_replica_dirs
export parse_llr, llr_dir_hdf5
export sort_by_central_energy_to_hdf5
include("llr_plots.jl")
export a_vs_central_action, full_trajectory_plot
include("clean_corrupted_replicas.jl")
export remove_non_matching_trajectories_in_replicas
include("histogram.jl")
export probability_density, plot_plaquette_histogram!
include("critical_histogram.jl")
include("double_gaussian_fit.jl")
include("free_energy.jl")
export plot_free_energies
include("provenance.jl")
export provenance, print_provenance_tex, print_provenance_csv, write_provenance_hdf5
include("supercooling.jl")
export supercooling
include("polyakov_loop.jl")
export is_fixed_a_measured, read_fixed_data
export logZ_fixed_a, polyakov_loop_fixed_a

end # module LLRParsing
