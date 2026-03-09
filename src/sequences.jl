"""
    ImageSequence{M}(f, images)

Lazily produces images of type `M`.
`images` may be an iterable collection of any data type `I` (e.g. `String` for filenames) so long as
`f(images)` returns a result of type `M`.
"""
struct ImageSequence{M<:AbstractMatrix,F,I} <: AbstractVector{M}
    f::F
    images::I
end

ImageSequence{M}(f::F, images::I) where {M,F,I} = ImageSequence{M,F,I}(f, vec(images))

Base.size(i::ImageSequence) = size(i.images)

@propagate_inbounds function Base.getindex(i::ImageSequence{M}, index::Integer) where M
    return i.f(i.images[index])::M
end

#---Define these methods so that image data isn't needlessly generated when printing---------------#

function Base.show(io::IO, i::ImageSequence{M}) where M
    print(io, ImageSequence{M}, (i.f, i.images))
end

function Base.summary(io::IO, i::ImageSequence{M,<:Any,<:AbstractVector{<:AbstractString}}) where M
    print(io, "Sequence of ", length(i), " images (type ", M, ") from filenames")
end

function Base.show(io::IO, ::MIME"text/plain", i::ImageSequence{M}) where M
    summary(io, i)
    isempty(i.images) && return nothing
    println(io, ':')
    # Seriously just a hack
    buf = IOBuffer()
    try
        show(IOContext(buf, :compact => true), MIME("text/plain"), i.images)
        seekstart(buf)
        readline(buf)
        write(io, buf)
    finally
        close(buf)
    end
end
