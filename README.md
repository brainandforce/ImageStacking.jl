# ImageStacking.jl

[![Build Status][ci-status-img]][ci-status-url]
[![Aqua.jl][aqua-img]][aqua-url]

This Julia package provides tools for performing astronomical image stacking.
The goal is to support commonly used image stacking methods that use robust means and efficiently apply them to large images with parallelization.
It also aims to provide other useful statistics that are not commonly used in astrophotography, such as the standard error of the mean, in hopes that they may be useful in other image processing tasks such as denoising.

Currently supported image stacking methods include:
  * mean and median
  * sigma clipping, MAD clipping, and Winsorized sigma clipping

Normalization (additive, multiplicative, and scaling) is also supported with IKSS measures as well as the median and MAD.
Stacking outputs the mean/median as well as the standard deviation and rejection counts in a stack, which can be used to estimate the standard error of the mean.

Future goals include:
  * Generalized extreme Studentized deviate test (GESDT) rejection
  * support for weighted means
  * support for images registered with [Drizzle (variable-pixel linear reconstruction)][drizzle-paper]

Image registration, plate solving, calculation of weights for weighted means, and other processing tasks are outside the scope of this package.
The goal is to build an ecosystem of small interoperating packages that are easily composed.

This package is intended to be generic with respect to image formats, though users will want to provide linear data in the vast majority of cases.
The package contains an extension for improved FITS support with [FITSIO.jl].

## Acknowledgements

This work would not be possible without the prior work of the Siril team.
[Siril] is a free astronomical image processing tool with support for both astrophotographic and scientific workflows, including photometry.
Special thanks to Cyril Richard [lock042] for answering all my questions about Siril in the Astrobiscuit Discord server.

Further acknowledgements go out to members of the Astrobiscuit and AURIC Discord servers, particularly [RuBoLo][RuBoLo-github] and [August Broe (xBrowi)][Browi-github].

[repo-url]:         https://github.com/brainandforce/ImageStacking.jl
[docs-stable-img]:  https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]:  https://brainandforce.github.io/ImageStacking.jl/stable
[docs-dev-img]:     https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:     https://brainandforce.github.io/ImageStacking.jl/dev
[ci-status-img]:    https://github.com/brainandforce/ImageStacking.jl/workflows/CI/badge.svg
[ci-status-url]:    https://github.com/brainandforce/ImageStacking.jl/actions
[aqua-img]:         https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]:         https://github.com/JuliaTesting/Aqua.jl
[codecov-img]:      https://codecov.io/gh/brainandforce/ImageStacking.jl/branch/main/graph/badge.svg
[codecov-url]:      https://codecov.io/gh/brainandforce/ImageStacking.jl/
[drizzle-paper]:    https://arxiv.org/abs/astro-ph/9808087
[FITSIO.jl]:        https://github.com/JuliaAstro/FITSIO.jl
[SIRIL]:            https://www.siril.org/
[lock042]:          https://github.com/lock042
[RuBoLo-github]:    https://github.com/RuBoLo
[Browi-github]:     https://github.com/xBrowi
[RuBoLo-instagram]: https://astrobin.com/users/
[Browi-astrobin]:   https://astrobin.com/users/AugustBroe