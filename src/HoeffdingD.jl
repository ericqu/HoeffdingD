module HoeffdingD

export hoeffdingd

using StatsBase
using LazyArrays ## for @~ to speed up Qi calculations

"""
    reject_independence(D, n, α)

    Compute the D-test of independence. 
    Compute the D-test of independence. 
        Given:
        a D (the calculated D measure),
        a n (the number of observations), 
        and a α (0 < α < 1) the desired level of significance.
    Return true to indicate rejection of H₀ hypothesis of independence.
    The implementation is based upon chapter 9 in the paper from Wassily Hoeffding (1948).

"""

function reject_independence(D, n, α)
    ρ = √((2 * (n^2 +5n -32)) / (9n * (n-1)*(n-3)*(n-4)*α)) 
    return D > ρ
end

"""
    hoeffdingd(a::AbstractVector{Union{Missing,T}}, b::AbstractVector{Union{Missing,T}}, pvalue=true, α=0.05) where {T <: Real}

Compute the Hoeffding measure of independence of the two continuous random variables. As defined by Hoeffding [A Non-Parametric Test of Independence](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-19/issue-4/A-Non-Parametric-Test-of-Independence/10.1214/aoms/1177730150.full) in chapter 5.

This is a wrapper to remove missing data by pair, and then call the core hoeffdingd to return the adequate value.

"""
function hoeffdingd(a::AbstractVector{Union{Missing,T}}, b::AbstractVector{Union{Missing,T}}) where {T <: Number}

    missa = findall(!ismissing, a)
    missb = findall(!ismissing, b)
    gindexes = findall(in(missa), missb)

    nma = convert(Vector{T}, view(a, view(missa, gindexes)))
    nmb = convert(Vector{T}, view(b, view(missb, gindexes)))

    return hoeffdingd(nma, nmb)
end

"""
    function hoeffdingd(a::AbstractVector{<:T}, b::AbstractVector{<:T}, α) where {T <: Number}

Compute the Hoeffding measure of independence of the two continuous random variables. As defined by Hoeffding [A Non-Parametric Test of Independence](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-19/issue-4/A-Non-Parametric-Test-of-Independence/10.1214/aoms/1177730150.full) in chapter 5.

Compute the Hoeffding measure of dependence (D) together with the result of the test of independence when α is presented.
Check that α is between 0 and 1. 
Return the results from the core hoeffdingd and the result of the independence test in a tuple.

"""
function hoeffdingd(a::AbstractVector{<:T}, b::AbstractVector{<:T}, α) where {T <: Number}
    if (α <= 0)  || (α >= 1)  
        throw(DomainError(α, "α must be: 0 < α < 1"))
    end
    D = hoeffdingd(a,b)
    N = size(a,1)
    return D , reject_independence(D, N, α)
end

"""
    function hoeffdingd(a::AbstractVector{<:T}, b::AbstractVector{<:T}) where {T <: Number}

Compute the Hoeffding measure of independence of the two continuous random variables. As defined by Hoeffding [A Non-Parametric Test of Independence](https://projecteuclid.org/journals/annals-of-mathematical-statistics/volume-19/issue-4/A-Non-Parametric-Test-of-Independence/10.1214/aoms/1177730150.full) in chapter 5.

Returns missing when there is not enough observations (less than 5).

The implementation relies on:
    StatsBase.tiedrank(a) to avoid duplication, 
    and LazyArrays (@~) to improve performance.

"""
function hoeffdingd(a::AbstractVector{<:T}, b::AbstractVector{<:T}) where {T <: Number}
    N = size(a, 1)
    if (N != size(b, 1))
        error("Both array should have the same number of (non missing) observations")
    end
    if N < 5 
        # error("Not enough observations ($(N) < 5) for Hoeffding's D")
        return missing
    end

    R = StatsBase.tiedrank(a)
    S = StatsBase.tiedrank(b)
    Q = Vector{Float64}(undef, N)

    @inbounds for i in 1:N
        Q[i] = 1 + sum(@~( (R .< R[i]) .& (S .< S[i]) ))
        Q[i] += 0.25 * (sum(@~((R .== R[i]) .& (S .== S[i]))) - 1)
        Q[i] += 0.5 * sum(@~((R .== R[i]) .& (S .< S[i]))) 
        Q[i] += 0.5 * sum(@~((R .< R[i]) .& (S .== S[i])))
    end

    D1 = sum((Q .- 1) .* (Q .- 2))
    D2 = sum((R .- 1) .* (R .- 2) .* (S .- 1) .* (S .- 2))
    D3 = sum((R .- 2) .* (S .- 2) .* (Q .- 1))
    D = 30 * ((N - 2) * (N - 3) * D1 + D2 - 2 * (N - 2) * D3) / (N * (N - 1) * (N - 2) * (N - 3) * (N - 4))

    return D
end

"""
    hoeffdingd(m::AbstractMatrix{<:T}, α) where {T <: Number}

Compute the Hoeffding measure of independence and the test of independence for each pair of columns in the matrix m.
Returns a matrix with the results.

"""
function hoeffdingd(m::AbstractMatrix{<:T}, α) where {T <: Number}

    cols = size(m)[2]
    results = Array{Tuple{Float64, Bool}}(undef, cols, cols)
    @inbounds for i in 1:cols
        results[i, i] = hoeffdingd(view(m, :, i), view(m, :, i), α)
        @inbounds for j in i + 1:cols
            results[i, j] = hoeffdingd(view(m, :, i), view(m, :, j), α)
            results[j, i] = results[i,j] 
        end
    end
    return results

end

"""
    hoeffdingd(m::AbstractMatrix{<:T}) where {T <: Number}

Compute the Hoeffding measure of independence for each pair of columns in the matrix m.
Returns a matrix with the results.

"""
function hoeffdingd(m::AbstractMatrix{<:T}) where {T <: Number}

    cols = size(m)[2]
    results = Array{Float64}(undef, cols, cols)
    @inbounds for i in 1:cols
        results[i, i] = hoeffdingd(view(m, :, i), view(m, :, i))
        @inbounds for j in i + 1:cols
            results[i, j] = hoeffdingd(view(m, :, i), view(m, :, j))
            results[j, i] = results[i,j] 
        end
    end
    return results

end


end 
