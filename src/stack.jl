#---Stacking functions-----------------------------------------------------------------------------#
"""
    ImageStacking.stack!(
        output::AbstractMatrix,
        block::ImageStacking.StackBlock,
        r::RejectionMethod,
        [coeffs::NormalizationCoefficients,]
        [op::Union{typeof(+),typeof(*)}]
    )

Stacks the data in `block` and writes the results to `output`.
The region indices of `block` must be contained within `output`; it will not be resized to
accomodate data which may be out of bounds.
Data in `block` may be modified depending on the choice of rejection method.

If normalization coefficients and a choice of normalization method (either `+` for additive or `*`
for multiplicative) are provided, then normalization will be applied to the images before they are
stacked.
"""
function stack!(
    output::AbstractMatrix,
    block::StackBlock,
    r::RejectionMethod,
    coeffs::NormalizationCoefficients,
    op::Union{typeof(+),typeof(*)}
)
    for (i,p) in zip(CartesianIndices(target_axes(block)), eachpixel(block))
        output[i] = pixel_stack!(p, r, coeffs, op)
    end
    return output
end

function stack!(
    output::AbstractMatrix,
    block::StackBlock,
    r::RejectionMethod
)
    for (i,p) in zip(CartesianIndices(target_axes(block)), eachpixel(block))
        output[i] = pixel_stack!(p, r)
    end
    return output
end

function Base.stack(
    blocks::AbstractVector{StackBlock{T}},
    image_dims::NTuple{2,Integer},
    r::RejectionMethod,
    coeffs::NormalizationCoefficients,
    op::Union{typeof(+),typeof(*)}
) where T
    output = Matrix{PixelStats{T}}(undef, image_dims)
    Threads.@threads for sb in blocks
        stack!(output, sb, r, coeffs, op)
    end
    return output
end

function Base.stack(
    images::ImageSequence{<:AbstractMatrix{T}},
    image_dims::NTuple{2,Integer},
    r::RejectionMethod,
    coeffs::NormalizationCoefficients,
    op::Union{typeof(+),typeof(*)}
) where T
    blocks = create_blocks(images, image_dims)
    return stack(blocks, image_dims, r, coeffs, op)
end

function Base.stack(
    images::ImageSequence{<:AbstractMatrix{T}},
    image_dims::NTuple{2,Integer},
    r::RejectionMethod,
    e::NormalizationEstimator,
    op::Union{typeof(+),typeof(*)}
) where T
    @info "Calculating normalization..."
    coeffs = get_normalization(median, e, images)
    @info "Stacking images..."
    blocks = create_blocks(images, image_dims)
    return stack(blocks, image_dims, r, coeffs, op)
end
