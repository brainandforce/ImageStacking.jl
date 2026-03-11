module ImageStacking

using Statistics
using SpecialFunctions

import Base: @propagate_inbounds
import Statistics: mean, median, median!

include("sequences.jl")
export ImageSequence
include("parallelization.jl")
export eachframe, eachpixel
include("statistics.jl")
include("normalization.jl")
export get_normalization
export MedianMADEstimator, IKSSEstimator, LocationOnly, MedianEstimator
include("rejection.jl")
export rejections, accepted
export NoRejection,SigmaClipping, MADClipping, WinsorizedSigmaClipping
include("stack.jl")

#---Stack images with statistics-------------------------------------------------------------------#

export RejectionMethod

end