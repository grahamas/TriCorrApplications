abstract type AbstractEEG end
get_signal(::AbstractEEG) = error("unimplemented")
get_times(::AbstractEEG) = error("unimplemented")
get_channel_names(::AbstractEEG) = error("unimplemented")

abstract type AbstractTimeseriesEEG <: AbstractEEG end
struct TimeseriesEEGv1{T,DATA<:AbstractMatrix{T},IND<:AbstractMatrix{T}} <: AbstractTimeseriesEEG
    times::Vector{T}
    data::DATA#Matrix{T}
    channel_names::Vector{String}
    indicators::IND#Matrix{T}
    indicator_names::Vector{String}
end
function get_signal(eeg::AbstractTimeseriesEEG)
    eeg.data
end
get_times(eeg::AbstractTimeseriesEEG) = eeg.times
get_channel_names(eeg::AbstractTimeseriesEEG) = eeg.channel_names

function plot_eeg_traces(eeg::AbstractEEG; labels=get_channel_names(eeg), std_max=nothing, downsample_factor=1, layout=:row)
    arr = NamedDimsArray{(:channel, :time)}(get_signal(eeg))
    arr = if std_max !== nothing
        arr = copy(arr)
        for i_channel ∈ 1:size(arr,1)
            timeseries = arr[i_channel, :]
            std_i = std(dropmissing(timeseries))
            arr[i_channel, abs.(timeseries) .> (std_max * std_i)] .= NaN
        end
        arr
    else
        arr
    end
    time = get_times(eeg)[begin:downsample_factor:end]
    arr = arr[:,begin:downsample_factor:end]
    @show size(arr)
    tbl = if labels === nothing
        Tables.table(arr', header=["$i" for i ∈ axes(arr, :channel)])
    else
        Tables.table(arr', header=labels)
    end
    df = DataFrame(tbl)
    df.time = time
    stacked_df = stack(df, axes(arr, :channel))
    rename!(stacked_df, :value => :signal, :variable => :channel)
    aog_data = data(stacked_df)
    plt = if layout == :row
        aog_data * mapping(:time, :signal, row=:channel) * visual(Lines)
    elseif layout == :layout
        aog_data * mapping(:time, :signal, layout=:channel) * visual(Lines)
    else
        error("bad layout kwarg")
    end
    return plt
end

function draw_eeg_traces!(fig, times, eeg; title=nothing, kwargs...)
    layout = GridLayout()
    signals = get_signal(eeg)
    labels = get_channel_names(eeg)
    for i_sig ∈ 1:size(signals, 1)
        label = (labels)[i_sig]
        layout[i_sig,1] = ax = Axis(fig)
        plot_contribution!(ax, eeg, times, signals[i_sig,:])
        layout[i_sig,2] = Label(fig, label, tellwidth=true, tellheight=false, rotation=-pi/2)
    end
    if title !== nothing
        layout[0,:] = Label(fig, title, tellwidth=false)
    end
    layout
end

function draw_eeg_traces!(fig, eeg; kwargs...)
    times = get_times(eeg)
    draw_eeg_traces!(fig, times, eeg; kwargs...)
end

function draw_eeg_traces(eeg; resolution, kwargs...)
    fig = Figure(resolution=resolution)
    layout = draw_eeg_traces!(fig, eeg; kwargs...) 
    fig[1,1] = layout
    fig
end

function plot_eeg_hists(maybe_arr::AbstractArray; labels=nothing)
    arr = NamedDimsArray{(:channel, :time)}(maybe_arr)
    tbl = if labels === nothing
        Tables.table(arr', header=["$i" for i ∈ axes(arr, :channel)])
    else
        Tables.table(arr', header=labels)
    end
    df = DataFrame(tbl)
    df.time = 1:nrow(df)
    stacked_df = stack(df, axes(arr, :channel))
    rename!(stacked_df, :value => :signal, :variable => :channel)
    aog_data = data(stacked_df)
    plt = aog_data * mapping(:signal, row=:channel) * AlgebraOfGraphics.histogram(bins=200)
    return plt
end

function draw_eeg_hists(arr::AbstractEEG; title=nothing, kwargs...)
    plt = plot_eeg_hists(arr; kwargs...)
    axis = (height=40, width=400)
    facet = (; linkyaxes = :none, linkxaxes=:none, grid=false, spines=false)
    fg = draw(plt; axis, facet)
    for axis_entry in fg.grid
        ax = axis_entry.axis
        tightlimits!(ax)
        #hidedecorations!(ax)
    end
    if title !== nothing
        Label(fg.figure[0,:], title, tellwidth=false)
    end
    fg
end

function plot_eeg_hists(eeg::AbstractEEG; kwargs...)
    plot_eeg_hists(eeg.signals; labels=eeg.labels, kwargs...)
end

using Makie

# ANY SCRIPT USING THIS MUST IMPORT A MAKIE
using DataFrames, AlgebraOfGraphics


plot_contributions(vec::AbstractVector; kwargs...) = plot_contributions(DataFrame([vec], :auto); kwargs...)
plot_contributions(arr::AbstractMatrix; kwargs...) = plot_contributions(DataFrame(arr', :auto); kwargs...)

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")



function plot_contributions!(fig, eeg::AbstractEEG, times::AbstractVector, data::AbstractArray; title=nothing, get_label, default_label_fns = Dict(
    "tricorr" => offset_motif_numeral,
    "aEEG" => (i -> eeg.labels[i])
), consensus_plot=true, other_params...)
    if get_label isa String
        get_label = default_label_fns[get_label]
    end
    n_rows = size(data, 1)
    layout = GridLayout()
    ax_plt_pairs = map(1:n_rows) do i_row
        layout[i_row, 1] = ax = Axis(fig, xlabel="", xgridvisible=false, ygridvisible=false)
        plt = plot_contribution!(ax, eeg, times, data[i_row,:]; other_params...)
        if !isnothing(get_label)
            signal_name=get_label(i_row)
            layout[i_row, 2] = Label(fig, signal_name, tellheight=false, tellwidth=true, rotation=-pi/2)
        end
        (ax, plt)
    end
    linkaxes!(first.(ax_plt_pairs)...)
    for (ax,plt) ∈ ax_plt_pairs[1:end-1]
        hidedecorations!(ax)
    end
    if consensus_plot
        layout[end+1,1] = cons_ax = Axis(fig)
        plot_reviewer_consensus!(cons_ax, eeg)
        linkxaxes!(first.(ax_plt_pairs)..., cons_ax)
    end
    if title !== nothing
        layout[0,:] = Label(fig, title, tellwidth=false)
    end
    return layout
end

function plot_contributions(eeg::AbstractEEG, times::AbstractVector, data::AbstractArray; resolution, kwargs...)
    fig = Figure(resolution = resolution, font = noto_sans)
    fig[1,1] = plot_contributions!(fig, eeg, times, data; kwargs...)
    return fig
end

function plot_contribution(eeg::AbstractEEG, args...; title=nothing, resolution=(800,600), kwargs...)
    fig = Figure(resolution=resolution)
    fig[1,1] = ax = Axis(fig)
    fig[2,1] = cons_ax = Axis(fig)
    l = plot_contribution!(ax, eeg, args...; kwargs...)
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    plot_reviewer_consensus!(cons_ax, eeg)
    linkxaxes!(ax, cons_ax)
    rowsize!(fig.layout, 3, Auto(0.3))
    return (fig, ax, l)
end

function plot_contribution!(ax::Axis, eeg, times, data; epoch_s=nothing)
    l = lines!(ax, times, data, color=:grey11, linecolor=:grey11)
    tightlimits!(ax); hidespines!(ax)
    hidedecorations!(ax, ticklabels=false)
    seizure_annotations, artifact_annotations = if isnothing(epoch_s)
        (eeg.seizure_annotations, eeg.artifact_annotations)
    else
        (
            discretize_and_merge_bounds(eeg.seizure_annotations,epoch_s),
            discretize_and_merge_bounds(eeg.artifact_annotations,epoch_s)
        )
    end
    seizure_onsets = [on for (on,off) in seizure_annotations if on < off]
    seizure_offsets = [off for (on,off) in seizure_annotations if on < off]
    artifact_onsets = [on for (on,off) in artifact_annotations if on < off]
    artifact_offsets = [off for (on,off) in artifact_annotations if on < off]
    if !isempty(artifact_onsets)
        vspan!(ax, artifact_onsets, artifact_offsets, color=(:red, 0.2))
    end
    if !isempty(seizure_onsets)
        vspan!(ax, seizure_onsets, seizure_offsets, color=(:blue, 0.2))
    end
    l
end

plot_reviewer_consensus!(ax, eeg::AbstractEEG) = error("unimplemented")
