function normalize_01!(arr)
    arr .-= minimum(arr)
    arr ./= maximum(arr)
    return arr
end

function seqclass_contributions_channels_white(snippet, boundary, λ_max...; n_bootstraps, bootstraps_step)
    f = fit(Whitening, snippet, dims=2)
    snippet .= MultivariateStats.transform(f, snippet)
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max...; n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step)
end

function calc_contributions_timeseries_channels_white(data, times, λ_n, λ_t; 
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
        contributions[:,i_window] = seqclass_contributions_channels_white(extended_window, boundary, λ_n, λ_t;
            n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
        )
        #ProgressMeter.next!(p)
    end
    return (contributions, times[window_starts])
end

function seqclass_contributions_channel_znorm(snippet, boundary, λ_max...; n_bootstraps, bootstraps_step)
    for ch ∈ axes(snippet, 1)
        snippet[ch,:] .-= mean(snippet[ch,:])
        snippet[ch,:] ./= std(snippet[ch,:])
    end
    normalize_01!(snippet)
    bootstrap_normed_sequence_classes(snippet, boundary, λ_max...; n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step)
end

function calc_contributions_timeseries_channel_znorm(data, times, λ_n, λ_t; 
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
        contributions[:,i_window] = seqclass_contributions_channel_znorm(extended_window, boundary, λ_n, λ_t;
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
    t0, t1 = boundary.t_bounds
    # normalize_01!(@view snippet[:,1:t0-1])
    # normalize_01!(@view snippet[:,t0:t1])
    # normalize_01!(@view snippet[:,t1+1:end])
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
    extended_window = Matrix{Float64}(undef, size(data,1), window_t + 2λ_t)
    @show typeof(extended_window) size(extended_window)
    map(1:10) do i_window#length(window_starts)) do i_window
        window_start = window_starts[i_window]
        window_end = window_start + window_t - 1
        copyto!(extended_window, @view data[:,window_start-λ_t:window_end+λ_t])
        ret_ext_win = copy(extended_window)
        boundary = PeriodicExtended(
            (λ_t+1,λ_t+window_t)
        )
        contributions[:,i_window] = seqclass_contributions_snippet_01(extended_window, boundary, λ_n, λ_t;
            n_bootstraps=n_bootstraps, bootstraps_step=bootstraps_step
        )
        ret_ext_win
        #ProgressMeter.next!(p)
    end
    #return (contributions, times[window_starts])
end