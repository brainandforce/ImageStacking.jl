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
If the data is stored in the most significant bits, we'll need to divide our ADU values by conversion factor, ``f``, to get the correct counts.
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
The mean is simply the sum of all data points in a sample ``X``, ``x_{k}``, divided by the number of points ``N``:
```math
\begin{aligned}
\operatorname{mean}\left(X\right) = \mu_{X} & = \sum_{k = 1}^{N} \frac{x_{k}}{N}
\end{aligned}
```

Along with the mean ``\mu_{X}``, we can also calculate the sample variance ``\sigma^{2}_{X} ``:
```math
\begin{aligned}
\operatorname{var}\left(X\right) = \sigma^{2}_{X} & = \sum_{k = 1}^{N} \frac{\left(x_{k} - \mu\right)^2}{N - 1}
\end{aligned}
```
While an unbiased estimator for the sample standard deviation has a relatively complicated form, we can take the square root of the sample variance to get a biased (but useful) estimate of the sample standard deviation.

Aside from the sample variance, we can also calculate the variance of the mean.
The sample variance and sample standard deviation tell us how spread out the individual measurement values are, but the variance of the mean tells us how accurate the estimate of the mean is.
[Bienaymé's identity](https://en.wikipedia.org/wiki/Bienaym%C3%A9%27s_identity) allows us to calculate the variance of the sum of variables from their individual variances:
```math
\begin{aligned}
\operatorname{var}\left(\sum_{k = 1}^{N} x_{k}\right) & = \sum_{k = 1}^{N} \operatorname{var}\left(x_{k}\right) + 2 \sum_{i = 1}^{N} \sum_{j = i + 1}^{N} \operatorname{cov}\left(x_{i}, x_{j}\right)
\end{aligned}
```
where ``\operatorname{cov}\left(x_{i}, x_{j}\right)`` is the covariance function, which is assumed to be zero if all measurements are independent.
Since the mean is just the sum of all ``x_{k}`` divided by the number of samples ``N``, the variance of the mean is simply the variance divided by ``N``:
```math
\begin{aligned}
\operatorname{var}\left(\mu_{x}\right) = \frac{\sigma^{2}_{X}}{N}
\end{aligned}
```
From this, we can estimate the standard deviation of the mean, which is more commonly called the *standard error of the mean*:
```math
\begin{aligned}
\operatorname{SE}\left(\mu_{x}\right) \approx \sqrt{\frac{\sigma^{2}_{X}}{N}} = \frac{\sigma_{X}}{\sqrt{N}}
\end{aligned}
```
The standard error of the mean decreases with the square root of the number of samples, so to halve the noise in a given image, the frame count must be quadrupled.

The problem with mean stacking is that the mean is not robust, with a breakdown point of 0.
Given a dataset of any size, it is possible to change the mean to an arbitrary value by altering the value of only a single datum.
Outliers are commonplace in astronomical imaging: aside from bad pixels and satellite trails, random cosmic ray strikes can often be seen on individual frames.
For this reason, mean stacking is generally not recommended unless the goal is to visualize artifacts.

### Weighted mean

It is also possible to define a weighted mean, given a set of weights ``w_{k}`` associated with data points ``x_{k}``:
```math
\begin{aligned}
\operatorname{mean}\left(X\right) = \mu_{X} & = \frac{\sum_{k = 1}^{N} w_{k} x_{k}}{\sum_{k = 1}^{N} w_k}
\end{aligned}
```

In the context of astronomical stacking, weights often correspond to the reliability or confidence in an observation.
If the variance ``\sigma^{2}_{k}`` associated with a particular data point ``x_{k}`` is known, ``w_{k}`` is simply the inverse of the variance.

## Median stacking

The median ``m_{X}`` of a dataset ``X`` with size ``N`` is a value for which exactly half of the dataset is greater than or equal to ``m_{X}`` and exactly half is less than or equal to ``m_{X}``:
```math
\begin{aligned}
m_{x} = \operatorname{median}\left(X\right) & = \frac{1}{2} \left(x_{\lceil \frac{N}{2} \rceil} + x_{\lfloor \frac{N + 1}{2} \rfloor}\right)
\end{aligned}
```
For odd ``N``, the median is the middle value when all data in ``X`` is sorted, since ``\lceil \frac{N}{2} \rceil`` and ``\lfloor \frac{N + 1}{2} \rfloor`` are identical.
For even ``N``, it is the average of the two middle values of sorted ``X``.

Unlike the mean, the median has a breakdown point of ``\frac{1}{2}``, which means that it can still yield a reasonable estimate of central tendency even if up to half of the dataset is arbitrarily corrupted.
This makes it is robust to artifacts that commonly appear in astronomical images.

Analogous to how the variance is the mean of the squared differences of each datum from the mean, the *median absolute deviation* is the median of the absolute differences of each datum from the median:
```math
\begin{aligned}
\operatorname{mad}\left(X\right) & = \operatorname{median}\left(\left| x_{k} - \operatorname{median}\left(X\right)\right|\right)
\end{aligned}
```

It is a bit less effective than the mean at reducing noise: the standard error in the median is approximately
```math
\operatorname{SE}\left(m_{X}\right) \approx \sqrt{\frac{2}{\tau}} = \frac{\sigma_{X}}{\sqrt{N}}
```
assuming that ``X`` is normally distributed.

The problem with median stacking is that it can only provide a maximum of 1 bit of extra precision in our estimate of the photon rate.
When working with a few images, this is not a problem, but as we increase our count of images, we cannot resolve differences in photon rate that are smaller than half of the camera's minimum ADU value.
This leads to posterization, especially in areas with faint detail, and it can cause banding artifacts in areas where the signal does not vary much.
As a result, median stacking is generally not recommended, even for small datasets, as there are alternative stacking methods which provide results of similar quality.

## Rejection stacking

The most useful image stacking methods combine the precision-increasing and noise-reducing properties of the mean with the robustness of the median.

### Trimmed mean

Given trimming fractions ``t_{\text{low}}`` and ``t_{\text{high}}`` in the interval ``[0, \frac{1}{2})`` (often chosen to be the same for a symmetric distribution), the *trimmed mean* (or *truncated mean*) is the mean of all elements of the dataset except for the ``t_{\text{low}} N `` smallest and ``t_{\text{high}} N`` largest elements.
The ordinary mean is recovered when ``t_{\text{low}} = t_{\text{high}} = 0``, and the median is recovered in the limit where ``t_{\text{low}} = t_{\text{high}} = \frac{1}{2}``.

In many cases, the trimming fractions do not divide the number of elements in the dataset evenly.
To allow the trimmed mean to be smoothly interpolated, the trimmed mean can be thought of as a weighted mean, with accepted points having weight 1 and discarded points having weight 1.
The smallest ``\lfloor t_{\text{low}} N \rfloor`` and largest ``\lfloor t_{\text{high}} N \rfloor`` elements will have weight zero.
The smallest element not discarded will have weight ``t_{\text{low}} N - \lfloor t_{\text{low}} N \rfloor``, and the largest element not discarded will have weight ``t_{\text{high}} N - \lfloor t_{\text{high}} N \rfloor``.

The standard error of the trimmed mean is the Winsorized standard deviation divided by the product of the number of data points and the untrimmed fraction of data.

While a trimmed mean is suitable for small datasets, it often discards good data.
For larger datasets, it's worth using methods that try to identify outliers rather than simply assuming a certain fraction of data is bad.

### Winsorized mean

Similar to the trimmed mean, the *Winsorized mean* uses Winsorizing fractions ``W_{\text{low}}`` and ``W_{\text{high}}``.
Instead of discarding these values, as in the trimmed mean, the values are replaced with their next nearest neighbors before the average is taken.

!!! todo
    How do we handle the cases where the Winsorizing fractions are not divisible by the number of data points?
    How can we calculate the standard error of a Winsorized mean?

### Sigma clipping

Given trimming parameters ``k_{\text{low}}`` and ``k_{\text{high}}``, this method calculates the median `m` and standard deviation ``\sigma`` of the data, then rejects any data below ``m - k_{\text{low}} \sigma`` and above ``m + k_{\text{high}} \sigma``.
This process is repeated on the trimmed data until no more points are trimmed.

### MAD clipping

This method is almost identical to sigma clipping, but uses the median absolute deviation rather than the standard deviation to determine whether to reject points.

!!! warning
    The MAD multiples used in this method are **not** scaled to match the standard deviation of a normal distribution.

### Winsorized sigma clipping

This method extends sigma clipping by adding an initial censorship step, limiting data to a fraction of the standard deviation from the median and recalculating the standard deviation until no pixels are rejected or the change in the standard deviation falls under a given threshold.
This new standard deviation is then used to perform sigma clipping.
As with sigma clipping, the process is repeated until no more pixels are rejected.

In other implementations of Winsorized sigma clipping, the censorship thresholds ``W_{\text{low}}`` and ``W_{\text{high}}`` are 1.5 standard deviations above and below the median, but they are adjustable in this implementation.
During each censorship iteration, the Winsorized standard deviation is scaled by a factor of 
```math
2 - \frac{1}{2} \left( \operatorname{erf}\left(\frac{W_{\text{low}}}{\sqrt{2}}\right) + \operatorname{erf}\left(\frac{W_{\text{high}}}{\sqrt{2}}\right) \right)
```
which is approximately `1.13361440253771617` when ``W_{\text{low}} = W_{\text{high}} = 1.5``.

!!! note
    Although this method is known as "Winsorized sigma clipping", the initial step is *not* Winsorization, as it censors the data based on standard deviations instead of using the low and high fractions.
    This name is currently used in PixInsight and SIRIL for this method, and it has been used here for consistency, but a different name may be used in the future.

### Generalized extreme Studentized deviate test (GESDT)

!!! todo
    Implement this method...
