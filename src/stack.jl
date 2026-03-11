function stack!(
    output::AbstractMatrix,
    block::StackBlock,
    r::RejectionMethod,
    coeffs::NormalizationCoefficients,
    op::Union{typeof(+),typeof(*)}
)
    for (i,p) in zip(CartesianIndices(block.region), eachpixel(block))
        output[i] = pixel_stack!(p, r, coeffs, op)
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
