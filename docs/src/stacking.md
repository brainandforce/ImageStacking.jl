# Stacking

## Theory

Even with a theoretically perfect sensor, it is impossible to acquire an image that is free of noise, because light is quantized.
The arrival of photons on a photosite follows a [Poisson distribution](https://en.wikipedia.org/wiki/Poisson_distribution) with probability mass function
```math
\begin{aligned}
\operatorname{Pois}(\lambda, k) & = \frac{\lambda^{k} e^{-\lambda}}{k!}
\end{aligned}
```
where ``\lambda`` is the rate parameter (expected number of photons per unit of time), ``k`` is the number of events, and the value of ``\operatorname{Pois}(\lambda, k)`` is the probability of measuring ``k`` events over the unit of time.
The mean (``\mu``) of ``\operatorname{Pois}(\lambda)``is ``\lambda``, as is its variance (``\sigma^{2}``), so its standard deviation (``\sigma``) is ``\sqrt{\lambda}``.
With a large enough value of ``\lambda``, we can treat a Poisson distribution as roughly equal to a Gaussian distribution with the same mean and square root.

When we take an image with a modern CMOS sensor, we obtain ADU values that correspond to the photon count, though indirectly.
Before taking an image, the sensor fills each of the CMOS wells with a pre-determined number of extra photons to prevent any data from clipping to zero.
This introduces an ADU bias or offset, ``B``, is set to prevent any pixels from having negative values, and the sensor operates at a certain gain setting.
While you may have provided your camera with a gain value (or ISO for more ordinary digital cameras), the number does not usually correspond to the conversion factor from photons to ADU, ``C_{\gamma}``.
However, it is possible to calculate your sensor's photon conversion factor at any gain setting by taking the ratio of the variance to the mean for a sequence of flat frames.

One more potential complication is that sensors may not use every bit of the data type used to store ADU values.
Many popular cameras for astrophotography are 14-bit cameras, so the ADU value is represented as a `UInt16`.
Depending on the software used, the camera may use the least significant or most significant bits to record the data.
If the data is stored in a e'll need to divide our ADU values by conversion factor, ``f``, to get the correct counts.
For a 14-bit sensor using the most significant bits of a `UInt16`, that factor is ``2^{16-14} = 4``.

With these numbers known, we can now estimate the number of photons collected in a single exposure, ``n_{\gamma}``, with a known ADU value ``A``:
```math
\begin{aligned}
n_{\gamma} & = \frac{A - B}{C_{\gamma} \times f}
\end{aligned}
```
This estimate is imperfect due to the presence of other noise sources, collectively known as read noise.
This includes thermal noise, random telegraph noise, quantization error (particularly when ``C_{\gamma} < 1``), and variation in the behavior of individual pixels across a sensor.

## Mean stacking

Although ``\sigma`` increases with the square root of ``\lambda``, the ratio ``\frac{\lambda}{\sigma}`` decreases with the square root of ``\lambda``.
This tells us that if we can collect more photons, our estimate of ``\lambda`` will likely be more accurate, as the potential spread in its values is smaller.
While this could theoretically be accomplished by taking a single exposure of extreme length, we can also sum or average many exposures to achieve the same goal.

Along with the mean ``\mu``, we can also calculate the sample variance ``s^{2}``:
```math
\begin{aligned}
s^{2} & = \sum_{k = 1}^{N} = \frac{\left(x_{k} - \mu\right)^2}{N - 1}
\end{aligned}
```
We can obtain the standard deviation

The problem with mean stacking is that the mean has a breakdown point of 0, so it is extremely sensitive to outliers.
Outliers are commonplace in astronomical imaging: aside from bad pixels and satellite trails, random cosmic ray strikes can often be seen on individual frames.

## Median stacking

Unlike the mean, the median has a breakdown point of ``\frac{1}{2}``, and it will get rid of any artifacts that are only present in a few images.

The problem with median stacking is that it can only provide a maximum of 1 bit of extra precision in our estimate of the photon rate.
When working with a few images, this is not a problem, but as we increase our count of images, we cannot resolve differences in photon rate that are smaller than half of the camera's minimum ADU value.
This leads to posterization, especially in areas with faint detail, and it prevents median stacking from bringing out faint detail.

## Rejection stacking

We can think about the medi

### Winsorized sigma clipping

This
