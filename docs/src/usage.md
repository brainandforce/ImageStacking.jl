# Usage

## Getting started

### Setup

This package is currently not available from the Julia registry, so you'll want to add this repo to your environment from the Git url:
```
(@v1.X) pkg> add https://github.com/esitohi/ImageStacking.jl.git
```

Additionally, you'll probably want to install [FITSIO.jl](https://github.com/JuliaAstro/FITSIO.jl) or some other package capable of reading and writing the format of the image files you wish to stack.

When using this package, we recommend ensuring that Julia is started with more than one thread, as the stacking operations take advantage of all available threads.
This may be manually reduced if desired.

### Stacking some data

To use ImageStacking.jl, you must have image data which has already been registered (or does not require registration, such as flat frames).
Image registration is beyond the scope of this package, but it may be accomplished with other software.
For this example, we assume you have data that has been calibrated and registered in Siril, so that it is associated with a sequence file (`*.seq`).
The function below will take a sequence file as input and return a `Vector{String}` of all the filenames of the sequence.
```julia
function filenames_from_seq(seqfile::AbstractString)
    dir = first(splitdir(seqfile))
    prefix = open(seqfile) do io
        readline(io)
        readline(io)
        readuntil(io, '\'')
        readuntil(io, '\'')
    end
    r = Regex(prefix * "\\d+\\.fits")
    return filter!(s -> !isnothing(match(r, s)), readdir(dir, join=true))
end
```
From the list of files, we can construct an `ImageSequence`.
By default, Siril works in 32-bit floating point, so change the `Float32` element if needed.
```julia-repl
julia> images = ImageSequence{Matrix{Float32}}(FITSIO.fitsread, filenames_from_seq("r_pp_light_.seq"))
```
In the future, this package may take a hard dependency on [AstronomicalImageSequences.jl](https://github.com/esitohi/AstronomicalImageSequences.jl), which should simplify the import of data.

Now that we have an `ImageSequence`, we can stack it.
For now, we will need to supply the dimensions of the final stacked image.
If all frames are the same size, this can simply be acquired from the size of the first frame:
```julia-repl
julia> image_dims = size(first(images))
```
Then we can supply all our arguments to `stack`: the image set `images` and the stacked image dimensions `image_dims`, along with our choice of stacking and normalization methods.
In most cases, we want to use a rejection stacking method, such as `WinsorizedSigmaClipping`.
We specify the rejection bounds for the sigma clipping step with two arguments

In Siril, the most common normalization estimator is the IKSS estimator, but location and scale estimates based on the median and MAD can also be used for faster normalization.
Both of these normalization estimators are available as `IKSSEstimator` and `MedianMADEstimator`.

When stacking light frames, we want to apply additive normalization with scaling, so our final argument should be `+`.
For flat frames, we only want to apply multiplicative normalization without scaling, so our final argument should be `*`, and we should wrap our estimator type in `LocationOnly` to omit the scaling correction.
For lights, the stacking call should look like this:
```julia-repl
julia> result = stack(images, image_dims, WinsorizedSigmaClipping(3,3), IKSSEstimator(), +)
```
The result contains much more information than the pixel values: the `PixelStats` data type also includes dispersion data (standard deviation), the number of pixels in the sample (which may be nonzero due to dithering around the edges), and the number of low and high rejections.
To extract the stacked image, use a list comprehension:
```julia-repl
julia> [x.scale for x in result]
```

!!! todo
    We are seriously lacking convenience functions here...

!!! todo
    Test if this works with frames that have been drizzled.
