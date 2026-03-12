module ImageStackingFITSIOExt

using ImageStacking
using FITSIO

import ImageStacking: ImageSequence, StackBlock

function ImageStacking.create_blocks(
    ::Type{T},
    filenames,
    image_dims::NTuple{2,Integer},
    memory_fraction::Real = 0.5,
    nthreads = Threads.nthreads(:default),
    free_memory = Sys.free_memory()
) where T<:Real
    memory_needed = sizeof(T) * prod(image_dims) * (length(filenames) + 1)
    working_memory = round(Int, free_memory * memory_fraction, RoundDown)
    if memory_needed < working_memory
        # Split the image up into chunks equal to the memory size
        (chunksize, remainder) = divrem(last(image_dims), nthreads, RoundNearest)
        regions = [
            (1:first(image_dims), ((n-1)*chunksize + 1):(n*chunksize + (n == nthreads) * remainder))
            for n in 1:nthreads
        ]
        blocks = [StackBlock{T}(undef, r, length(filenames)) for r in regions]
        for (i, fn) in enumerate(filenames)
            FITS(fn, "r", extendedparser=false) do f
                for b in blocks
                    read!(f[1], view(b.buffer, :, :, i), b.region...)
                end
            end
        end
    else
        ImageStacking._insufficient_memory_error(memory_needed, working_memory, memory_fraction)
    end
    return blocks
end

end