
noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")

function plot_djia_stocks(djia=get_stocks_df(djia_tickers))
    plt = data(djia) * mapping(:date, :close, color=:ticker) * visual(Lines)
end

function plot_contributions(dates, contributions; title=nothing, resolution, events=nothing, get_label=offset_motif_numeral)
    n_motifs = size(contributions, 1)
    fig = Figure(resolution=resolution, font=noto_sans)
    map(1:n_motifs) do i_row
        ax = Axis(fig[i_row, 1], xlabel="", xgridvisible=false, ygridvisible=false)
        plt = plot_contribution!(ax, dates, contributions[i_row,:]; 
            events=events
        )
        motif=get_label(i_row)
        Label(fig[i_row, 2], motif, tellheight=false, tellwidth=true, rotation=-pi/2)
        (ax, plt)
    end
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    return fig
end

function plot_contribution!(ax, times::AbstractVector{<:Date}, contribution::AbstractVector; events=nothing)
    @warn "Plotting dates unsupported; defaulting to ints."
    l = lines!(ax, 1:length(times), contribution)
    tightlimits!(ax); hidespines!(ax)
    if events !== nothing
        onset_dxs = [findlast(Ref(on) .>= times) for (on, off) ∈ events if times[1] <= on < off < times[end]]
        offset_dxs = [findlast(Ref(off) .>= times) for (on, off) ∈ events if  times[1] <= on < off < times[end]]
        filter!(!isnothing, onset_dxs); filter!(!isnothing, offset_dxs)
        onsets = collect(1:length(times))[onset_dxs]
        offsets = collect(1:length(times))[offset_dxs]
        if !isempty(onsets)
            vspan!(ax, onsets, offsets, color=(:grey, 0.2))
        end
    end
    l
end

function plot_contribution(args...; title=nothing, resolution=(800,600), kwargs...)
    fig = Figure(resolution=resolution)
    ax = Axis(fig[1,1])
    l = plot_contribution!(ax, args...; kwargs...)
    if title !== nothing
        Label(fig[0,:], title, tellwidth=false)
    end
    return (fig, ax, l)
end

# https://www.rosettacode.org/wiki/Roman_numerals/Encode#Julia
function roman_encode(n::Integer)
    if n < 1 || n > 4999 throw(DomainError(n)) end
 
    DR = [["I", "X", "C", "M"] ["V", "L", "D", "MMM"]]
    rnum = ""
    for (omag, d) in enumerate(digits(n))
        if d == 0
            omr = ""
        elseif d <  4
            omr = DR[omag, 1] ^ d
        elseif d == 4
            omr = DR[omag, 1] * DR[omag, 2]
        elseif d == 5
            omr = DR[omag, 2]
        elseif d <  9
            omr = DR[omag, 2] * DR[omag, 1] ^ (d - 5)
        else
            omr = DR[omag, 1] * DR[omag + 1, 1]
        end
        rnum = omr * rnum
    end
    return rnum
end