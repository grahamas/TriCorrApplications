function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end

function seqclass_contributions_channel_white(snippet, boundary, λ_max...; n_bootstraps, bootstraps_step)
    for ch ∈ axes(snippet, 1)
        snippet[ch,:] .-= mean(snippet[ch,:])
        snippet[ch,:] ./= std(snippet[ch,:])
    end
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max...; n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step)
end

function calc_contributions_timeseries_channel_white(data, times, λ_n, λ_t; 
        window_t = 2λ_t+1, 
        n_bootstraps=2, bootstraps_step=2
    )
    n_days = size(data, 2)
    window_starts = (λ_t+1):window_t:(n_days-window_t-λ_t)
    contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Union{Float64,Missing}, 14, length(window_starts)))
    #p = ProgressMeter.Progress(length(window_starts))
    for i_window ∈ 1:length(window_starts)
        window_start = window_starts[i_window]
        window_end = window_start + window_t - 1
        extended_window = data[:,window_start-λ_t:window_end+λ_t] 
        boundary = PeriodicExtended(
            (λ_t+1,λ_t+window_t)
        )
        contributions[:,i_window] = seqclass_contributions_channel_white(extended_window, boundary, λ_n, λ_t;
            n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
        )
        #ProgressMeter.next!(p)
    end
    return (contributions, times[window_starts])
end

function seqclass_contributions_channel_01(snippet, boundary, λ_max...; n_bootstraps, bootstraps_step)
    for ch ∈ axes(snippet, 1)
        normalize_01!(snippet[ch,:])
    end
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max...; n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step)
end

function calc_contributions_timeseries_channel_01(data, times, λ_n, λ_t; 
        window_t = 2λ_t+1, 
        n_bootstraps=2, bootstraps_step=2
    )
    n_days = size(data, 2)
    window_starts = (λ_t+1):window_t:(n_days-window_t-λ_t)
    contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Union{Float64,Missing}, 14, length(window_starts)))
    #p = ProgressMeter.Progress(length(window_starts))
    for i_window ∈ 1:length(window_starts)
        window_start = window_starts[i_window]
        window_end = window_start + window_t - 1
        extended_window = data[:,window_start-λ_t:window_end+λ_t] 
        boundary = PeriodicExtended(
            (λ_t+1,λ_t+window_t)
        )
        contributions[:,i_window] = seqclass_contributions_channel_01(extended_window, boundary, λ_n, λ_t;
            n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
        )
        #ProgressMeter.next!(p)
    end
    return (contributions, times[window_starts])
end

function seqclass_contributions_snippet_01(snippet, boundary, λ_max...; n_bootstraps, bootstraps_step)
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max...; n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step)
end

function calc_contributions_timeseries_snippet_01(data, times, λ_n, λ_t; 
        window_t = 2λ_t+1, 
        n_bootstraps=2, bootstraps_step=2
    )
    n_days = size(data, 2)
    window_starts = (λ_t+1):window_t:(n_days-window_t-λ_t)
    contributions = NamedDimsArray{(:motif_class, :time)}(zeros(Union{Float64,Missing}, 14, length(window_starts)))
    #p = ProgressMeter.Progress(length(window_starts))
    for i_window ∈ 1:length(window_starts)
        window_start = window_starts[i_window]
        window_end = window_start + window_t - 1
        extended_window = data[:,window_start-λ_t:window_end+λ_t] 
        boundary = PeriodicExtended(
            (λ_t+1,λ_t+window_t)
        )
        contributions[:,i_window] = seqclass_contributions_snippet_01(extended_window, boundary, λ_n, λ_t;
            n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
        )
        #ProgressMeter.next!(p)
    end
    return (contributions, times[window_starts])
end