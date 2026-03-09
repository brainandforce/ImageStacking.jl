"""
    NormalizationCoefficients{T<:Real}

Stores the offset and scale data from a normalization method.
The median of the offset and scale are cached
"""
struct NormalizationCoefficients{T<:Real} <: AbstractVector{Tuple{T,T}}
    data::Vector{Tuple{T,T}}
    centers::Tuple{T,T}
end

function NormalizationCoefficients(f, locations::AbstractVector, dispersions::AbstractVector)
    centers = promote(f(locations), f(dispersions))
    T = promote_type(eltype(centers), eltype(locations), eltype(dispersions))
    return NormalizationCoefficients{T}(collect(T, zip(locations, dispersions)), centers)
end

function NormalizationCoefficients(
    index::Integer,
    locations::AbstractVector,
    dispersions::AbstractVector
)
    centers = promote(locations[index], dispersions[index])
    T = promote_type(eltype(locations), eltype(dispersions))
    return NormalizationCoefficients{T}(collect(T, zip(locations, dispersions)), centers)
end

function NormalizationCoefficients(f, locations::AbstractVector)
    center = f(locations)
    T = promote_type(typeof(center), eltype(locations))
    return NormalizationCoefficients{T}([(l, one(T)) for l in locations], (center, one(T)))
end

function NormalizationCoefficients(index::Integer, locations::AbstractVector{T}) where T
    center = locations[index]
    return NormalizationCoefficients{T}([(l, one(T)) for l in locations], (center, one(T)))
end

Base.size(c::NormalizationCoefficients) = size(c.data)

@propagate_inbounds function Base.getindex(c::NormalizationCoefficients, i::Int)
    return (first(c.data[i]), last(c.data[i]))
end

"""
    NormalizationEstimator{T}

Represents a method of calculating the location and dispersion of an image, with both results
represented as type `T`.
To apply a normalization estimator, call an instance of it as a function.
"""
abstract type NormalizationEstimator{T}
end

"""
    MedianMADEstimator{T} <: NormalizationEstimator{T}

Calculates location with the median and scale with the median absolute deviation (MAD).
"""
struct MedianMADEstimator{T} <: NormalizationEstimator{T}
end

function (::MedianMADEstimator)(pixels)
    p = sort!(vec(copy(pixels)))
    return (mediansorted(p), madsorted(p))
end

"""
    IKSSEstimator{T} <: NormalizationEstimator{T}

Iterative k-sigma estimator of location and scale.
"""
struct IKSSEstimator{T} <: NormalizationEstimator{T}
end

(::IKSSEstimator{T})(pixels; k = 4, tol = eps(Float32)) where T = convert.(T, ikss(pixels; k, tol))

struct LocationOnly{T,E<:NormalizationEstimator{T}} <: NormalizationEstimator{T}
end

function (::LocationOnly{E,T})(pixels; kwargs...) where {E,T}
    return (convert(T, first(E()(pixels; kwargs...))), one(T))
end

const MedianEstimator{T} = LocationOnly{MedianMADEstimator{T},T}

(::MedianEstimator{T})(pixels) where T = (convert(T, median(pixels)), one(T))

#---Get normalization coefficients for an image sequence-------------------------------------------#
"""
    get_normalization(f, e::NormalizationEstimator{T}, images) -> NormalizationCoefficents{T}
    get_normalization(i::Integer, e::NormalizationEstimator{T}, images)

Calculates normalization coefficients for a sequence of images.
The normalizations may be provided with reference to a central estimate of the location and
dispersion parameters through a function `f`, or to the index `i` of a reference image.

Because normalization can take a long time for each image depending on the method used, this
function is multithreaded and will use all available threads to compute the normalization of a
sequence.
"""
function get_normalization(f, e::NormalizationEstimator{T}, images) where T
    ld = Vector{Tuple{T,T}}(undef, length(images))
    Threads.@threads for i in eachindex(images)
        ld[i] = e(images[i])
    end
    return NormalizationCoefficients{T}(ld, (f(first.(ld)), f(last.(ld))))
end

function get_normalization(i::Integer, e::NormalizationEstimator{T}, images) where T
    ld = Vector{Tuple{T,T}}(undef, length(images))
    Threads.@threads for i in eachindex(images)
        ld[i] = e(images[i])
    end
    return NormalizationCoefficients{T}(ld, (first.(ld[i]), last(ld[i])))
end

#---Perform normalization on entire image sequences------------------------------------------------#
"""
    apply_normalization!(
        op::Union{typeof(+),typeof(*)},
        pixel!::AbstractVector,
        coeffs::ImageStacking.NormalizationCoefficents
    )

Applies normalization to a set of pixel values.
This may be additive or multiplicative, depending on whether `op` is `+` or `*`.
"""
function apply_normalization!(
    ::typeof(+),
    pixel!::AbstractVector,
    coeffs::NormalizationCoefficients
)
    (l0, d0) = coeffs.centers
    for i in eachindex(pixel!)
        (l, d) = coeffs[i]
        pixel![i] = muladd(pixel![i] - l, d0 / d, l0)
    end
    return pixel!
end

function apply_normalization!(
    ::typeof(*),
    pixel!::AbstractVector,
    ndata::AbstractVector{<:NormalizationCoefficients}
)
    (l0, d0) = coeffs.centers
    for i in eachindex(pixel!)
        (l, d) = coeffs[i]
        pixel![i] = pixel![i] * (l0 / l) * (d0 / d)
    end
    return pixel!
end
