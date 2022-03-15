module TriCorrApplications

using CSV, DataFrames, DataFramesMeta,
    Dates,
    CairoMakie, 
    AlgebraOfGraphics,
    ProgressMeter,
    NamedDims, Statistics, StatsBase, DSP

using TripleCorrelations

include("../data/djia_tickers.jl")
include("../data/recessions.jl")

include("load.jl")
export get_clean_djia, timeseries_from_df,
    djia_tickers, recessions_since_1971

include("plot.jl")
export plot_contributions, plot_contribution

include("contributions.jl")
export calc_contributions_timeseries_snippet_01,
	calc_contributions_timeseries_channel_01,
	calc_contributions_timeseries_channel_white

end
