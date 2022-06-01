using Base.Threads
using DrWatson
@quickactivate "TriCorrApplications"

using CairoMakie
using AlgebraOfGraphics
# using DSP, Statistics, StatsBase
ext = "png"
using Random, JLD2, Dates, NamedDims, Statistics, DataFrames, MAT

using TriCorrApplications, TripleCorrelations

let (λ_n, λ_t) = (8,10), #recalculate=true,
    n_bootstraps=0, bootstraps_step=2,
    start_date = Date("1981-05-01",dateformat"yyyy-mm-dd"),
    task_name = "djia_DEBUG_start_$(Dates.format(start_date,"yyyy_mm_dd"))_lags_$(λ_n)_$(λ_t)",
    task_time = Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS");

boundary = PeriodicExtended((11,31))
debug_rasters = matread(datadir("debug_rasters_djia.mat"))["debug_rasters"]
# map(debug_rasters[2:end]) do raster
#     TriCorrApplications.seqclass_contributions_snippet_01(raster, boundary, λ_n, λ_t;
#             n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
#         )
# end

# djia_df = get_clean_djia(start_date=start_date)
# (djia, trading_dates) = timeseries_from_df(djia_df)
# djia

end