module TriCorrApplications

using CSV, DataFrames, DataFramesMeta,
    Dates,
    CairoMakie, 
    AlgebraOfGraphics,
    ProgressMeter,
    NamedDims, Statistics, StatsBase, DSP,
    MultivariateStats,
    NamedDims, Dates,
    MAT


using TripleCorrelations

include("../data/djia_tickers.jl")
include("../data/recessions.jl")

include("load/stocks.jl")
export get_clean_djia, timeseries_from_df,
    djia_tickers, recessions_since_1971,
    load_cattan_alpha_subject

include("util.jl")
export discretize_bounds, merge_bounds, discretize_and_merge_bounds
include("plot.jl")
export plot_contributions, plot_contribution,
    plot_reviewer_consensus!

include("contributions.jl")
export calc_contributions_timeseries_snippet_01,
	calc_contributions_timeseries_channel_01,
	calc_contributions_timeseries_channel_znorm,
    calc_contributions_timeseries_channels_white

include("eeg.jl")
export AbstractEEG, get_signal, 
    get_times, get_channel_names,
    plot_contributions, plot_eeg_traces, draw_eeg_traces,
    draw_eeg_traces!, plot_contributions!

end
