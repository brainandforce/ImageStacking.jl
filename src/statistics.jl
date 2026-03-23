"""
    NORMAL_STD_SQRTBWMV_RATIO

The ratio betwen the standard deviation and the square root of the biweight midvariance of a
standard normal distribution.
"""
const NORMAL_STD_SQRTBWMV_RATIO = 0.9909048195090278

#---Median of a sorted array-----------------------------------------------------------------------#

@inline function mediansorted(A::AbstractVector)
    inds = ((firstindex(A) + lastindex(A)) >>> 1) .+ (0:iseven(length(A)))
    return middle(view(A, inds))
end

#---Median absolute deviation----------------------------------------------------------------------#
"""
    madm(itr, m)

Compute the median absolute deviation (MAD) of a collection `itr`, with the known median `m`.
"""
madm(itr, m) = @inline median(abs(x - m) for x in itr)

# TODO: sorted version?
@inline function madsorted(A::AbstractVector)
    m = mediansorted(A)
    deviations = copy(A)
    for i in eachindex(deviations)
        deviations[i] = abs(deviations[i] - m)
    end
    # deviations = [abs(x - m) for x in A]
    T = eltype(deviations)
    # Deviations array is sorted in reverse up to this index
    s = (firstindex(A) + lastindex(A)) >>> 1
    # to ensure X is smaller than Y, index from 1:s-1
    # if A has even length, A[s] = A[s+1], so the subarrays are still sorted
    # if A has odd length, A[s] = 0 (absolute deviations must be positive)
    # since we made our own vector through comprehension, the indexing does not need to be generic
    X = @view deviations[reverse(1:s-1)]
    Y = @view deviations[s:end]
    # faster to avoid evaluating these lengths in the loop
    lenX, lenY = length(X), length(Y)
    half = (lenX + lenY) >>> 1
    l, r = (1, s-1)
    # TODO: refactor this because it could be much simpler!
    # The view X does not actually need to be reversed
    while true
        i = (l + r) >>> 1
        # defined so that 
        j = half - i
        Xl = (i >= 1) ? (@inbounds X[i]) : typemin(T)
        Xr = (i <= lenX) ? (@inbounds X[i+1]) : typemax(T)
        Yl = (j >= 1) ? (@inbounds Y[j]) : typemin(T)
        Yr = (j <= lenY) ? (@inbounds Y[j+1]) : typemax(T)
        if (Xl <= Yr) && (Yl <= Xr)
            tmp = min(Xr, Yr)
            return isodd(lenX + lenY) ? tmp : middle(max(Xl, Yl), tmp)
        elseif Xl > Yr
            r = i - 1
        elseif Yl > Xr
            l = i + 1
        end
    end
end

"""
    mad(itr)

Compute the median absolute deviation (MAD) of a collection `itr`.
"""
mad(itr) = @inline madm(itr, median(itr))

#---Biweight midvariance---------------------------------------------------------------------------#
"""
    bwmvm(itr, m)

Compute the biweight midvariance (BWMV) of the collection `itr` given a known median `m`.
"""
@inline function bwmvm(itr, m)
    d = madm(itr, m)
    # Use the eltype of the MAD for type stability reasons
    numerator = denominator = zero(eltype(d))
    for x in itr
        u = abs(x - m) / (9 * d)
        if u < 1
            w = (1 - u^2)
            numerator += ((x - m) * w^2)^2
            denominator += w*(1 - 5*u^2)
        end
    end
    return length(itr) * numerator/(denominator^2)
end

"""
    bwmvmsorted(itr, m)

Compute the biweight midvariance (BWMV) of a sorted collection `itr` more efficiently.
"""
@inline function bwmvsorted(itr)
    m = mediansorted(itr)
    d = madsorted(itr)
    numerator = denominator = zero(eltype(d))
    for x in itr
        u = abs(x - m) / (9 * d)    # 9 * sqrt(2) * erfinv(1/2) times the stdev if normal
        if u < 1
            w = (1 - u^2)
            numerator += ((x - m) * w^2)^2
            denominator += w*(1 - 5*u^2)
        end
    end
    return length(itr) * numerator/(denominator^2)
end

"""
    bwmv(itr)

Compute the biweight midvariance of the collection `itr`.
"""
bwmv(itr) = bwmvm(itr, median(itr))

#---Iterative k-sigma estimation-------------------------------------------------------------------#
"""
    ikss!(A, k = 4, tol = eps(Float32)) -> Tuple{Real, Real, Int, Int}

Performs iterative k-sigma estimation of location and scale (IKSS) on an array `A`, which is
modified in-place.
This estimator is used for image normalization and scaling.

The parameters of the output are:
* location
* scale
* count of pixels not rejected
* number of iterations performed
"""
function ikss!(A::AbstractVector; k = 4, tol = eps(Float32))
    P = sort!(A)
    i,j = (firstindex(P), lastindex(P))
    # For type stability
    s0 = middle(one(eltype(P)))
    for it in Iterators.countfrom(1)
        if isempty(i:j)
            return (zero(s0), zero(s0))
        end
        v = @view P[i:j]
        # the array is already sorted, so a full median is wasteful
        m = mediansorted(v)
        s = sqrt(bwmvsorted(v))
        if (iszero(s) || abs(s0 - s) / s < tol)
            return (m, convert(eltype(P), s * NORMAL_STD_SQRTBWMV_RATIO))
        end
        # Throw out the samples above and below
        i = findfirst(>=(m - k*s), P)
        j = findlast(<=(m + k*s), P)
        s0 = s
    end
end

ikss!(A::AbstractArray; k = 4, tol = eps(Float32)) = ikss!(vec(A); k, tol)

"""
    ikss(itr, k = 4, tol = eps(Float32)) -> Tuple{Real, Real, Int, Int}

Performs iterative k-sigma estimation of location and scale (IKSS).
This estimator is used for image normalization and scaling.

The parameters of the output are:
* location
* scale
* count of pixels not rejected
* number of iterations performed
"""
ikss(itr; k = 4, tol = eps(Float32)) = ikss!(collect(itr); k, tol)

#---Means with rejection---------------------------------------------------------------------------#
"""
    trimmed_mean_weights(T::Type, sz::Integer, lo::Real, hi::Real)

Calculates weights needed for the trimmed mean with trimming fractions `lo` and `hi`, handling cases
where the size of the dataset `sz` multiplied by either `lo` or `hi` is not an integer.
"""
function trimmed_mean_weights(T::Type, sz::Integer, lo::Real, hi::Real)
    weights = ones(T, sz)
    # Get weights for the lowest values
    extent_lo = lo * sz
    n_discard_lo = floor(Int, extent_lo)
    weights[begin:n_discard_lo] .= 0
    weights[n_discard_lo + 1] .= extent_lo - n_discard_lo
    # Get weights for the highest values
    extent_hi = hi * sz
    n_discard_hi = floor(Int, extent_hi)
    weights[end - n_discard_hi] = extent_hi - n_discard_hi
    weights[end - n_discard_hi + 1:end] .= 0
    return weights
end

# TODO: maybe create a lazy data structure for this?

function trimmed_mean!(data::AbstractVector{T}, lo::Real, hi::Real) where T
    sort!(data) # TODO: this should probably be a partialsort!
    weights = trimmed_mean_weights(float(T), length(data), lo, hi)
    return sum(x*w for (x, w) in zip(weights, data))
end

trimmed_mean!(data::AbstractVector, lohi::Real) = trimmed_mean!(data, lohi, lohi)

trimmed_mean(data::AbstractVector, lo::Real, hi::Real) = trimmed_mean!(copy(data), lo, hi)
trimmed_mean(data::AbstractVector, lohi::Real) = trimmed_mean!(copy(data), lohi, lohi)

# NOTE: the code above is probably not good to use in the general case
# We can just calculate the weight vector once and reuse it constantly

#=
"""
    winsorize_by_sigma!(A!::AbstractVector, l, h)
    winsorize_by_sigma!(A!::AbstractVector, k)

Winsorizes the data in `A!` which is `l` times the standard deviation of `A!` below the mean and `h`
times the standard deviation of `A!` above the mean, or `k` standard deviations away from the mean
in either direction if only one argument is given.
"""
function winsorize_by_sigma!(A!::AbstractVector, l, h)
    m = mean(A!)
    s = stdm(A!, m)
    w_lo = m - l*s
    w_hi = m + h*s
    for i in eachindex(A!)
        A![i] = min(w_hi, max(w_lo, A![i]))
    end
    return A!
end

winsorize_by_sigma!(A!::AbstractVector, lh) = winsorize_by_sigma(A!, lh, lh)
=#

#---Pixel statistics-------------------------------------------------------------------------------#
"""
    PixelStats{T}

Represents statistical information for a pixel in a stack:
  * The scale of the data (usually the mean)
  * The dispersion of the data (usually the standard deviation, but it is the median absolute 
    deviation for `R<:MedianRejection`)
  * The number of pixels in the sample
  * The number of pixels rejected for being too low and too high
"""
struct PixelStats{T}
    scale::T
    dispersion::T
    count::Int
    reject_lo::Int
    reject_hi::Int
end

function PixelStats(
    scale::T,
    dispersion::T,
    count::Integer,
    reject_lo::Integer,
    reject_hi::Integer
) where T
    return PixelStats{T}(rejection, scale, dispersion, count, reject_lo, reject_hi)
end

function PixelStats(scale, dispersion, count::Integer, reject_lo::Integer, reject_hi::Integer)
    T = promote_type(typeof(mean), typeof(stdev))
    return PixelStats{T}(scale, dispersion, count, reject_lo, reject_hi)
end

"""
    rejections(p::PixelStats) -> Tuple{Int,Int}

Returns the number of pixels rejected for falling below or above the rejection thresholds.
"""
rejections(p::PixelStats) = (p.reject_lo, p.reject_hi)

"""
    accepted(p::PixelStats) -> Int

Returns the number of pixels accepted in the stack.
This is equal to `p.count - sum(rejections(p))`.
"""
accepted(p::PixelStats) = p.count - sum(rejections(p))

"""
    std_err_mean(p::PixelStats)

Estimates the standard error in the mean from the standard deviation and number of accepted pixels
in the stack.
"""
std_err_mean(p::PixelStats) = p.dispersion / sqrt(accepted(p))