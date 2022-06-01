# Load financial data (src: https://stooq.com/db/h/)

function get_stock_df(ticker::String, walk_itr=walkdir(joinpath(@__DIR__, "..", "data", "exp_raw", "d_us_txt")))
    ticker = lowercase(ticker)
    ticker, country_ext = splitext(ticker)
    ticker = if country_ext == ""
        "$(ticker).us"
    else
        "$(ticker)$(country_ext)"
    end
    stock_file = reduce(vcat, map(walk_itr) do (root, dirs, files)
        joinpath.(Ref(root), filter(files) do file
            file == "$(ticker).txt"
        end)
    end) |> only
    df = CSV.read(stock_file, DataFrame)
    select(df, "<TICKER>" => :ticker, "<DATE>" => (x -> Date.(string.(x), Ref(dateformat"yyyymmdd"))) => :date, "<CLOSE>" => :close)
end

function get_stocks_df(tickers::Vector{String})
    df = mapreduce(get_stock_df, vcat, tickers)
end

function get_stocks_mx(tickers::Vector{String})
    (Matrix ∘ pre_mx_df ∘ get_stocks_df)(tickers)
end

function pre_mx_df(df::DataFrame)
    sort!(stack(df, :date, :ticker, :close), :date)
end

count_not_missing(x) = count(.!ismissing.(x))
function get_clean_djia(; start_date, max_n_days_missing=100)
    df = get_stocks_df(djia_tickers)
    @rsubset!(df, :date >= start_date)
    sufficient_trading_days = @chain df begin
        @by(:ticker, :n_trading_days = count_not_missing(:close))
        @subset(:n_trading_days .>= (maximum(:n_trading_days) - max_n_days_missing))
    end
    @rsubset!(df, :ticker ∈ sufficient_trading_days.ticker)
    df = unstack(df, :date, :ticker, :close)
    sort!(df, :date)
    dropmissing!(df)
    df
end

function timeseries_from_df(df::DataFrame)
    dates = df.date |> Vector
    timeseries = df[:, Not(:date)] |> Matrix
    (timeseries', dates)
end