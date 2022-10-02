
function offset_motif_numeral(n::Integer)
    roman_encode_zero(n-1)
end

# https://www.rosettacode.org/wiki/Roman_numerals/Encode#Julia
function roman_encode_zero(n::Integer)
    if n == 0 return "0" end
    if n < 0 || n > 4999 throw(DomainError(n)) end
 
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

function step_below(num::T, step, num_0) where T
    floor(Int, (num - num_0) / step) * step
end
function step_above(num::T, step, num_0) where T
    ceil(Int, (num - num_0) / step) * step
end

function discretize_bounds(bounds::ARR, step, bound_0=0; min_bound, max_bound) where ARR
    # FIXME warning does not keep bounds within any constraints
    ARR(map(bounds) do (start, stop)
        (max(min_bound,step_below(start, step, bound_0)), min(max_bound,step_above(stop, step, bound_0)))
    end)
end

function merge_bounds(bounds1::AbstractVector{Tup}, bounds2::AbstractVector{Tup}) where {T,Tup<:Tuple{T,T}}
    if isempty(bounds1) && isempty(bounds2)
        return Tup[]
    elseif isempty(bounds1)
        return bounds2
    elseif isempty(bounds2)
        return bounds1
    end
    bounds = sort([bounds1..., bounds2...])
    new_bounds = Tup[]

    current_start, current_stop = first(bounds)
    for (start, stop) in bounds[2:end]
        if start <= current_stop
            current_stop = max(stop, current_stop)
        else
            push!(new_bounds, (current_start, current_stop))
            current_start = start; current_stop = stop
        end
    end
    push!(new_bounds, (current_start, current_stop))

    return new_bounds
end

function discretize_and_merge_bounds(bounds, step; min_bound, max_bound)
    posterized_bounds = discretize_bounds(bounds, step; min_bound=min_bound, max_bound=max_bound)
    return merge_bounds(posterized_bounds, posterized_bounds)
end