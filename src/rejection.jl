#---Stacking methods-------------------------------------------------------------------------------#
"""
    StackingMethod

Supertype for all methods used to stack images.

Subtypes of `StackingMethod` should include all parameters needed to perform the operation,
such as rejection thresholds, that do not depend on the images being stacked.
"""
abstract type StackingMethod
end

"""
    MeanStacking <: StackingMethod

Produces the mean without any outlier rejection.
This method is not recommended for real astronomical data, since outliers are common enough to
seriously affect the final result.
"""
struct MeanStacking <: StackingMethod
end

"""
    MedianStacking <: StackingMethod

Produces the median, which is essentially a fully clipped mean.
This method is not recommended for stacking any more than a few frames, since the median cannot
significantly increase the precision of the data.
This leads to posterization, especially in darker areas of astronomical images.
"""
struct MedianStacking <: StackingMethod
end

"""
    AbstractSigmaClipping{T} <: StackingMethod

Supertype for methods based on sigma clipping, which iteratively reject outliers lying a given 
number of standard deviations from the median until none are rejected.
"""
abstract type AbstractSigmaClipping{T<:Real} <: StackingMethod
end

reject_low(r::AbstractSigmaClipping) = r.lo
reject_high(r::AbstractSigmaClipping) = r.hi

"""
    SigmaClipping{T<:Real} <: AbstractSigmaClipping{T<:Real}

Ordinary sigma clipping based on the median and standard deviation of a dataset.
"""
struct SigmaClipping{T<:Real} <: AbstractSigmaClipping{T}
    lo::T
    hi::T
end

SigmaClipping(reject::T) where T = SigmaClipping{T}(reject, reject)
SigmaClipping{T}(reject::Real) where T = SigmaClipping{T}(reject, reject)

"""
    MADClipping{T<:Real} <: AbstractSigmaClipping{T<:Real}

Modified sigma clipping based on the median and median absolute deviation (MAD) of a dataset.
"""
struct MADClipping{T<:Real} <: AbstractSigmaClipping{T}
    lo::T
    hi::T
end

MADClipping(reject::T) where T = MADClipping{T}(reject, reject)
MADClipping{T}(reject::Real) where T = MADClipping{T}(reject, reject)

"""
    WinsorizedSigmaClipping{T<:Real} <: AbstractSigmaClipping{T<:Real}

Winsorized sigma clipping, a modified version of sigma clipping which first performs winsorization
on the data before sigma clipping.
This procedure is iterated until no samples are rejected.
"""
struct WinsorizedSigmaClipping{T<:Real} <: AbstractSigmaClipping{T}
    lo::T
    hi::T
end

WinsorizedSigmaClipping{T}(reject::Real) where T = WinsorizedSigmaClipping{T}(reject, reject)
WinsorizedSigmaClipping(reject::T) where T = WinsorizedSigmaClipping{T}(reject, reject)

function WinsorizedSigmaClipping(lo::Real, hi::Real)
    T = promote_type(typeof(lo), typeof(hi))
    return WinsorizedSigmaClipping{T}(lo, hi)
end

"""
    GeneralizedESD{T<:Real} <: StackingMethod

Performs the [generalized extreme Studentized deviate test (GESDT)][Rosner1983] to reject outliers.

[Rosner1983]: https://doi.org/10.1080/00401706.1983.10487848
"""
struct GeneralizedESD{T<:Real} <: StackingMethod
    fraction::T
end

#---Means with rejection---------------------------------------------------------------------------#
"""
    ImageStacking.pixel_stack!(A!::AbstractVector, s::StackingMethod) -> PixelStats{T}

Stacks the pixel data `A!` using the stacking method `s`, potentially modifying `A!` in place.
This returns a [`PixelStats`](@ref) object, containing the scale and dispersion of the data, along
with the number of sampled and rejected pixels.
Note that this function may not modify `A!` depending on the method used (such as calculating a 
mean without rejection).

This function should only be called if it is both possible and acceptable to modify the input.
Otherwise, use [`pixel_stack`](@ref).
"""
function pixel_stack!(A!::AbstractVector, ::MeanStacking)
    m = mean(A!)
    s = stdm(A!, m)
    return PixelStats(m, s, length(A!), 0, 0)
end

function pixel_stack!(A!::AbstractVector, ::MedianStacking)
    sort!(A!)
    m = mediansorted(A!)
    s = madsorted(A!)
    l = length(A!)
    c = 1 + iseven(l)   # only one or two pixels are used, really
    r = (l - c) >>> 1   # we're rejecting everything other than the central data points
    return PixelStats(m, s, c, r, r)
end

"""
    ImageStacking.sigma_clip!(A!::AbstractVector, s_lo, s_hi; by = std, sorted::Bool = false)
        -> Tuple{Int, Int}

Performs sigma clipping on `A!` after sorting it in place, returning the number of pixels rejected
for lying below and above the thresholds relative to the median.
By default, the standard deviation is used, but this can be changed to a more robust measure of 
statistical dispersion by changing the `by` parameter to a different function.
Ideally, the function used to calculate the dispersion at each iteration should be able to take
advantage of the sorting of the data (e.g. `ImageStacking.madsorted` is a better choice than
`ImageStacking.mad` for calculating the median absolute deviation).
"""
function sigma_clip!(dispersion, A!::AbstractVector, s_lo, s_hi; sorted::Bool = false)
    l = length(A!)
    sorted || sort!(A!)
    lo = 0  # number of low pixels rejected
    hi = 0  # number of high pixels rejected
    # iterations = 0
    num_rejections = 1  # number of rejections per cycle
    while (num_rejections > 0) && (l - (lo + hi) > 3)
        # iterations += 1
        data = @view A![begin+lo:end-hi]
        m = mediansorted(data)
        deviation = dispersion(data)
        # Calculate the thresholds for rejection
        r_lo = m - s_lo * deviation
        r_hi = m + s_hi * deviation
        num_rejections = 0
        for i in 0:(length(data) >>> 1)
            is_lo = data[begin+i] < r_lo
            is_hi = data[end-i] > r_hi
            num_rejections += (is_lo + is_hi)
            (!is_lo && !is_hi) && break
            lo += is_lo
            hi += is_hi
        end
    end
    return (lo, hi)
end

function sigma_clip!(dispersion, A!::AbstractVector, s_lo, s_hi, m, sigma; sorted::Bool = false)
    l = length(A!)
    sorted || sort!(A!)
    lo = 0  # number of low pixels rejected
    hi = 0  # number of high pixels rejected
    # iterations = 0
    num_rejections = 1  # number of rejections per cycle
    while l - (lo + hi) > 3
        # iterations += 1
        data = @view A![begin+lo:end-hi]
        # Calculate the thresholds for rejection
        r_lo = m - s_lo * sigma
        r_hi = m + s_hi * sigma
        num_rejections = 0
        for i in 0:(length(data) >>> 1)
            is_lo = data[begin+i] < r_lo
            is_hi = data[end-i] > r_hi
            num_rejections += (is_lo + is_hi)
            (!is_lo && !is_hi) && break
            lo += is_lo
            hi += is_hi
        end
        num_rejections > 0 || break
        m = mediansorted(data)
        sigma = dispersion(data)
    end
    return (lo, hi)
end

function sigma_clip!(A!::AbstractVector, s_lo, s_hi; by = std, sorted::Bool = false)
    return sigma_clip!(by, A!, s_lo, s_hi; sorted)
end

function sigma_clip!(A!::AbstractVector, s_lo, s_hi, m, sigma; by = std, sorted::Bool = false)
    return sigma_clip!(by, A!, s_lo, s_hi, m, sigma; sorted)
end

function pixel_stack!(A!::AbstractVector, r::SigmaClipping)
    (lo, hi) = sigma_clip!(A!, reject_low(r), reject_high(r))
    data = @view A![begin+lo:end-hi]
    m = mean(data)
    s = stdm(data, m)
    return PixelStats(m, s, length(A!), lo, hi)
end

function pixel_stack!(A!::AbstractVector, r::MADClipping)
    (lo, hi) = sigma_clip!(A!, reject_low(r), reject_high(r), by=madsorted)
    data = @view A![begin+lo:end-hi]
    m = mean(data)
    s = stdm(data, m)
    return PixelStats(m, s, length(A!), lo, hi)
end

function winsorize!(A!::AbstractVector, m, sigma, s_low, s_high; sorted::Bool = false)
    l = length(A!)
    sorted || sort!(A!)
    w_lo = m - sigma * s_low
    w_hi = m - sigma * s_high
    num_winsorized = 0
    for i in 0:(l >>> 1)
        is_lo = A![begin+i] < w_lo
        is_hi = A![end-i] > w_hi
        (!is_lo && !is_hi) && break
        num_winsorized += (is_lo + is_hi)
    end
    return num_winsorized
end

"""
    ImageStacking.winsorize_for_sigma_clip!(
        A!::AbstractVector;
        tol = 0.0005,
        s_low = 1.5,
        s_high = 1.5,
        sorted = false
    )

Performs iterative Winsorization of the input data with respect to its median and standard
deviation, sorting and modifying it in-place.
"""
function winsorize_for_sigma_clip!(
    A!::AbstractVector;
    tol = 0.0005,
    s_low = 3//2,
    s_high = 3//2,
    sorted::Bool = false
)
    l = length(A!)
    sorted || sort!(A!)
    m = mediansorted(A!)
    s = std(A!)
    s0 = zero(s)
    # This is where the factor of 1.134 comes in in other implementations
    std_factor = 2 - (erf(s_low / sqrt(2)) + erf(s_high / sqrt(2))) / 2
    while ((abs(s - s0) / s0) > tol) #=&& (num_rejections > 0) && (l - (lo + hi) > 3)=#
        num_winsorized = winsorize!(A!, m, s, s_low, s_high, sorted=true)
        num_winsorized > 0 || break
        m = mediansorted(A!)
        s0 = s
        s = std(A!) * std_factor
    end
    return (m, s)
end

function pixel_stack!(A!::AbstractVector, r::WinsorizedSigmaClipping)
    # Perform iterative winsorization and sigma clipping until no more changes occur
    sort!(A!)
    (lo, hi) = (0, 0)   # number of low/high pixels rejected
    inds = eachindex(A!)
    while length(inds) > 3
        data = @view A![inds]
        # Track the number of winsorized and sigma clipped samples
        (m, σ) = winsorize_for_sigma_clip!(data, sorted=true)
        (sl, sh) = sigma_clip!(data, reject_low(r), reject_high(r), m, σ, sorted=true)
        # Winsorized samples aren't considered rejected
        lo += sl
        hi += sh
        inds_old = inds
        inds = (first(inds) + sl):(last(inds) - sh)
        inds_old == inds && break
    end
    data = @view A![inds]
    m = mean(data)
    s = stdm(data, m)
    return PixelStats(m, s, length(A!), lo, hi)
end

"""
    pixel_stack!(
        A!::AbstractVector,
        s::StackingMethod,
        coeffs::NormalizationCoefficients,
        op::Union{typeof(+),typeof(*)}
    ) -> PixelStats{T}

Performs pixel stacking with normalization, given a set of normalization coefficients and either the
`+` or `*` operator depending on whether additive or multiplicative normalization is desired.
"""
function pixel_stack!(
    A!::AbstractVector,
    s::StackingMethod,
    coeffs::NormalizationCoefficients,
    op::Union{typeof(+),typeof(*)}
)
    return pixel_stack!(apply_normalization!(op, A!, coeffs), s)
end

"""
    pixel_stack([::Type{T}], itr, s::StackingMethod) -> PixelStats{R,T}

Stacks pixel data in an iterator `itr` with stacking method `s`.
This returns a [`PixelStats`](@ref) object, containing the scale and dispersion of the data, along
with the number of sampled and rejected pixels.

If it is both possible and non-disruptive to modify the pixel data in place, it may be preferable to
call [`ImageStacking.pixel_stack!`] instead, which does not allocate a new array.
"""
pixel_stack(::Type{T}, itr, s::StackingMethod) where T = pixel_stack!(collect(T, itr), s)
pixel_stack(itr, s::StackingMethod) = pixel_stack!(collect(itr), s)