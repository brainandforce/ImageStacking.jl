# Normalization

## Why normalize?

When calculating robust means, fluctuations in background illumination or transparency can either cause inadvertent rejection of good data or prevent the rejection of bad data.
If a frame is affected by a passing cirrus cloud or a nearby bright light turning on, the background light will increase the average pixel value throughout the frame, and it may be enough to cause pixels from the frame to rejected as outliers.
Alternatively, if one frame sequence is taken in a light polluted area and another is taken at a pristine dark site, the calculated standard deviation of each pixel may be inflated by the difference between the pixel values at the two sites rather than the expected variation from shot noise.

## Measures of location and scale

Normalization requires estimates of the location (average pixel value) and scale (dispersion of average pixel values).
These estimators need to be robust enough that stars do not significantly affect the estimate.
While their presence may not be a significant problem if images are taken with exactly identical framing, dithering may cause stars to move in and out of frame, which would significantly affect the estimate.

### Median/MAD

One simple robust way to estimate location and scale is to calculate the median and median absolute deviation (MAD) of each image.
While more statistically efficient methods are available, they are very computationally efficient, and may be suitable for use with enormous datasets.

### Biweight midvariance (BWMV)

The *biweight midvariance* (BWMV) is a robust estimator of scale that is more efficient than the MAD.

### IKSS estimator

The *iterative k-sigma estimator of location and scale* (IKSS) was available in PixInsight and is currently used in Siril to perform normalization.

## Normalization techniques

For light frames, we want to compensate for variations in background illumination and sky transparency, so additive normalization with scaling is recommended.
Even on a perfectly clear night, atmospheric absorption will decrease and increase as a target rises and sets, and the sky background can vary due to local light domes or even airglow at dark sites.

For flat frames, multiplicative normalization should be used.
Flat frames are intended for division, and should not contain any additive terms.
It is essential that a offset level, bias frame, or dark frame is subtracted from the raw flats before they are normalized and stacked.
Since flat frames are usually taken under very consistent conditions, scaling corrections are not needed.

For dark and bias frames, normalization should never be used, as the raw ADU values are critically important for correct calibration.

## Image weighting

Because shot noise is proportional to the square root of the measured signal, images with higher average illumination will have more intrinsic noise.
When stacking the images, we can weight them by the inverse of the location estimate to favor those expected to have less noise.
