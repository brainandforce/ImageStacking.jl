"""
    StackBlock{T<:Number} <: AbstractVector{Matrix{T<:Number}}

Contains all of the data needed to perform stacking.

`StackBlock` can only store numerical data.
Color data should be broken down into individual data channels before stacking.
"""
struct StackBlock{T<:Number} <: AbstractVector{Matrix{T}}
    # Region of each input image used to construct this block
    region::NTuple{2,UnitRange{Int}}
    # Used when calculating the mean, defaults to one(T) for all elements
    # weights::Vector{T}
    # Stores all of the image data
    buffer::Array{T,3}
end

function StackBlock{T}(
    ::UndefInitializer,
    region::NTuple{2,AbstractUnitRange{<:Integer}},
    image_count::Integer
) where T
    return StackBlock{T}(region, Array{T,3}(undef, length.(region)..., image_count))
end

"""
    StackBlock{T}(f, images, region::NTuple{2,R})

Generates a `StackBlock` from a collection of images.

The function `f` is used to transform the elements of `images` into data that can be stored in an
`AbstractMatrix`
"""
function StackBlock{T}(
    f,
    images,
    region::NTuple{2,R}
) where {T,R<:UnitRange}
    buffer = Array{T,3}(undef, size.(region)..., length(images))
    for (n,i) in enumerate(images)
        buffer[:,:,n] .= f(i)
    end
    return StackBlock{T}(region, buffer)
end

# Index the StackBlock like a vector of images
Base.size(sb::StackBlock) = tuple(size(sb.buffer, 3))
@propagate_inbounds Base.getindex(sb::StackBlock, i::Int) = sb.buffer[:,:,i]
@propagate_inbounds Base.setindex!(sb::StackBlock, x, i::Int) = (sb.buffer[:,:,i] = x)

"""
    target_axes(sb::StackBlock)

Returns the region of the output image that stacking of this block will write to.
"""
target_axes(sb::StackBlock) = sb.region

"""
    ImageStacking.eachframe(sb::StackBlock)

Creates slices of each image in a `StackBlock`.
"""
eachframe(sb::StackBlock) = eachslice(sb.buffer, dims = 3)

"""
    ImageStacking.eachpixel(sb::StackBlock)

Creates slices of each pixel in a `StackBlock`.
Operations that calculate pixel statistics should use this iterator.
"""
eachpixel(sb::StackBlock) = eachslice(sb.buffer, dims = (1,2))

function Base.summary(io::IO, sb::StackBlock)
    print(io, typeof(sb), " with ", length(sb), " images (pixel region ", target_axes(sb), ')')
end

#---Break down images into blocks for stacking-----------------------------------------------------#

function _insufficient_memory_error(memory_needed, working_memory, memory_fraction)
    @error LazyString(
        "Not enough memory is available to construct blocks for the whole image.",
        "\nTotal memory required: ", memory_needed >>> 20, " MiB",
        "\nTotal memory available: ", working_memory >>> 20, " MiB",
        "\nHave you tried increasing memory_fraction (currently ", memory_fraction, "?)"
    )
    throw(OutOfMemoryError())
end

"""
    create_blocks(
        T,
        image_size,
        image_count,
        [memory_fraction = 0.5],
        [nthreads = Threads.nthreads(:default),]
        [free_memory = Sys.free_memory]
    ) -> Vector{StackBlock{T}}

Creates a vector of uninitialized `StackBlock` objects that can hold the desired number of images of
the given size and data type.
The number of blocks depends on the number of threads desired, which can be determined automatically
or manually specified.
The `memory_fraction` parameter determines how much of the available memory
"""
function create_blocks(
    ::Type{T},
    image_dims::NTuple{2,Integer},
    image_count::Integer,
    memory_fraction::Real = 0.5,
    nthreads = Threads.nthreads(:default),
    free_memory = Sys.free_memory()
) where T<:Real
    # Pad memory requirements to account for the destination image and
    memory_needed = sizeof(T) * prod(image_dims) * (image_count + 2)
    working_memory = round(Int, free_memory * memory_fraction, RoundDown)
    if memory_needed < working_memory
        # Split the image up into chunks equal to the memory size
        (chunksize, remainder) = divrem(last(image_dims), nthreads, RoundNearest)
        regions = [
            (1:first(image_dims), ((n-1)*chunksize + 1):(n*chunksize + (n == nthreads) * remainder))
            for n in 1:nthreads
        ]
        return [StackBlock{T}(undef, r, image_count) for r in regions]
    else
        _insufficient_memory_error(memory_needed, working_memory, memory_fraction)
    end
end

function create_blocks(
    images::ImageSequence{<:AbstractMatrix{T}},
    image_dims::NTuple{2,Integer},
    memory_fraction::Real = 0.5,
    nthreads = Threads.nthreads(:default),
    free_memory = Sys.free_memory()
) where T<:Real
    image_count = length(images)
    memory_needed = sizeof(T) * prod(image_dims) * (image_count + 1)
    working_memory = round(Int, free_memory * memory_fraction, RoundDown)
    if memory_needed < working_memory
        # Split the image up into chunks equal to the memory size
        (chunksize, remainder) = divrem(last(image_dims), nthreads, RoundNearest)
        regions = [
            (1:first(image_dims), ((n-1)*chunksize + 1):(n*chunksize + (n == nthreads) * remainder))
            for n in 1:nthreads
        ]
        blocks = [StackBlock{T}(undef, r, length(images)) for r in regions]
        Threads.@threads for i in eachindex(images)
            image = images[i]
            for b in blocks
                b.buffer[:,:,i] .= view(image, target_axes(b)...)
            end
        end
        return blocks
    else
        _insufficient_memory_error(memory_needed, working_memory, memory_fraction)
    end
end
