

#10.1080/03610926.2018.1563164 and https://www.stat.cmu.edu/tr/tr776/tr776.pdf
# helper functions
Z(λ,ν) = exp(ν*λ^(1/ν))/((λ^((ν-1)/(2*ν)))*((2*π)^((ν-1)/2)) *sqrt(ν))
#com-poissson pmf
cppmf(loglambda,ν,y,logz) = y * loglambda - logz - ν * SpecialFunctions.loggamma(y+1)

struct CPM_Poisson <: DiscreteUnivariateDistribution
    λ
    ν
end
# this is a relatively slow but accurate way to draw from cpm-poisson.
# Shouldn't matter, as mostly we will be evaluationg the logpmf, not gereating random samples.
function Distributions.rand(rng::AbstractRNG, d::CPM_Poisson) 
    P = rand()
    y = -1
    p = 0
    logz = log(Z(d.λ,d.ν))
    loglambda = log(d.λ)
    while p < P
        y += 1
        p += exp(cppmf(loglambda,d.ν,y,logz))
    end
    return y
end


function Distributions.logpdf(d::CPM_Poisson, x::Real) 
    logz = log(Z(d.λ,d.ν))
    loglambda = log(d.λ)
    cppmf(loglambda,d.ν,x,logz)
end
