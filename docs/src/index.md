```@meta
CurrentModule = ImageStacking
```
# ImageStacking.jl

[ImageStacking.jl](https://github.com/brainandforce/ImageStacking.jl) provides utilities for stacking raw images.
This is an key operation in astronomical imaging, and there are many different ways to go about it depending on the imager's goals.

While there are many programs, free and proprietary, that can stack raw images, ImageStacking.jl has three major goals for providing unique functionality:
  * **Provide the finest possible control over stacking:** Some stacking methods use parameters that are not exposed to the user in popular software, and we want to allow users to alter whatever hidden parameters exist.
  * **Provide ways of obtaining more statistical information than just the mean:** of greatest interest is calculation of the standard errors, which may be useful for denoising.
  * **Thoroughly documenting all stacking methods provided in this package:** some of the more popular methods for stacking images are poorly documented, and this package can serve as both a reference implementation and an educational guide to how these methods work.

```@autodocs
Modules = [ImageStacking]
```