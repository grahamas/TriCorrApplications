using Base.Threads
using DrWatson
@quickactivate "TriCorrApplications"

using CairoMakie, AlgebraOfGraphics
# using DSP, Statistics, StatsBase
ext = "png"
using Random, JLD2, Dates, NamedDims, Statistics, DataFrames

using TriCorrApplications

let (λ_n, λ_t) = (8,10), recalculate=true,
    start_date = Date("1981-05-01",dateformat"yyyy-mm-dd"),
    task_name = "djia_seqcont_channels_white_start_$(Dates.format(start_date,"yyyy_mm_dd"))_lags_$(λ_n)_$(λ_t)",
    task_time = Dates.format(Dates.now(), "yyyy_mm_dd-HHMMSS");

plots_subdir = plotsdir("$(task_name)_$(task_time)_channels_white")
mkpath(plots_subdir)

contributions, contributions_dates = if recalculate == true
    djia_df = get_clean_djia(start_date=start_date)
    djia_fig = data(stack(djia_df, Not(:date))) * mapping(:date, :value, color=:variable) * visual(Lines) |> draw
    save(joinpath(plots_subdir, "djia_closing_prices.png"), djia_fig)
    (djia, trading_dates) = timeseries_from_df(djia_df)
    contributions, contributions_dates = calc_contributions_timeseries_channels_white(djia, trading_dates, λ_n, λ_t)
    save(datadir("exp_pro", "$(task_name)_$(task_time).jld2"), Dict("contributions" => contributions, "contributions_dates" => contributions_dates))
    (contributions, contributions_dates)
else
    all_saves = mapreduce(vcat, walkdir(datadir("exp_pro"))) do (root, dirs, files)
        basenames = filter(file -> occursin(task_name, file), files)
        joinpath.(Ref(root), basenames)
    end
    recent_save = sort!(all_saves)[end]
    @load recent_save contributions contributions_dates
    (contributions, contributions_dates)
end

conts_fig = plot_contributions(contributions_dates, contributions; title="Motif Contributions (DJIA)", resolution=(1000,1600), events=recessions_since_1971)
save(joinpath(plots_subdir, "djia_contributions.png"), conts_fig)

rms(xs) = sqrt(mean(xs .^ 2))
rms_timeseries = [rms(contributions[:,i_sec]) for i_sec ∈ 1:size(contributions,2)]
(rms_fig, ax, l) = plot_contribution(contributions_dates, rms_timeseries; events=recessions_since_1971, resolution=(1200, 500), title="DJIA RMS motif contributions (grey = recession)")
save(joinpath(plots_subdir, "djia_rms_contributions.png"), rms_fig)

end
