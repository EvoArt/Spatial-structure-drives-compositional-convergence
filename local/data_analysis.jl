#load packages
using CSV, DataFrames, DataFrameMacros,Turing,StatsBase,
 Colors,Chain,  SpecialFunctions,Random, Distances,nMDS, BetaDispersion
include("DiversityFunctions.jl")

include("cpmpoisson.jl")
n_iter = 2000 # number of MCMC samples to take per chain
n_keep = 1000 # number of MCMC samples to keep per chain
n_chn = 4 # MCMC chains per model

Random.seed!(123)

### Load data
# Load in and manipulate data
Gdata =CSV.read("denoised_q.csv", DataFrame)
n_gen = 17
X = Array(Gdata[:,2:n_gen]) # Array of genus counts
Y = X .> 0 # Presence/absence matrix
props = X ./ Gdata.depth # proportions
n_row,n_genera = size(X)
n_treat = length(unique(Gdata.treatment))
source_inds = [findfirst(x -> x == t,["ara","glu","leu"]) for t in Gdata.carbon];
transfer_inds = [findfirst(x -> x == t,["none","liquid","bead"]) for t in Gdata.transfer];
### Alpha diversity

### model definition
# We use a com-poisson distribution, since our counts are underdispersed
#10.1080/03610926.2018.1563164
# helper functions
#com-poissson pmf
# define model
# we use a hierachical model and parametize in terms of the mean and dispersion parameter
@model function hier_CP(Y,carbon_inds,σ,μ)

    mu ~ filldist(Normal(0,1),3)
    mu .*= σ
    mu .+= μ # carbon source means
    ν ~ truncated(Normal(0,5),lower = 0)

    if sum((mu .+((ν -1) ./(2 *ν))) .> 0) == 3
        λ = (mu .+((ν -1) ./(2 *ν))) .^(ν)
        Y .~  CPM_Poisson.(λ[carbon_inds],ν)
    else
        Turing.@addlogprob! -Inf
    end
end

@model function hier_CP2(d,e,f,g,h,i)
    sig ~ Exponential(5)
    muy ~ filldist(truncated(Normal(10,10),lower = 1),3)
    @submodel prefix="sub1" a = hier_CP(d,e,sig,muy[1])
    @submodel prefix="sub2" b = hier_CP(f,g,sig,muy[2])
    @submodel prefix="sub3" c = hier_CP(h,i,sig,muy[3])
end
mod2 = hier_CP2(Gdata.richness[Gdata.transfer .== "none"],source_inds[Gdata.transfer .== "none"],
                Gdata.richness[Gdata.transfer .== "liquid"],source_inds[Gdata.transfer .== "liquid"],
                Gdata.richness[Gdata.transfer .== "bead"],source_inds[Gdata.transfer .== "bead"])

chn = sample(mod2, NUTS(0.7),MCMCThreads(),n_iter,n_chn,drop_warmup = true)[n_iter-n_keep+1:end]     
CSV.write("richness_chn.csv",DataFrame(chn)) 

# evenness

Gdata.pielou = pielou.(eachrow(X))

@model function pielou_sep(x)
    
    p ~ filldist(truncated(Normal(0.5,1),0,1),3)
    dh ~ Exponential()
    d ~ filldist(Exponential(),3)
    c ~ arraydist([filldist(Normal(p[i],dh),3) for i in 1:3])
    for i in 1:3
        for j in 1:3
        x[i][j] .~ Normal(c[j,i],d[i])
        end
    end
end

sep_mod = pielou_sep([[Gdata.pielou[Gdata.treatment .== trt] for trt in ["ara_none","glu_none","leu_none"]],
[Gdata.pielou[Gdata.treatment .== trt] for trt in ["ara_liquid","glu_liquid","leu_liquid"]],
[Gdata.pielou[Gdata.treatment .== trt] for trt in ["ara_bead","glu_bead","leu_bead"]]])

sep_chn = sample(sep_mod, NUTS(0.7),MCMCThreads(),n_iter,n_chn,drop_warmup = true)[n_iter-n_keep+1:end]
CSV.write("pielou_chn.csv",DataFrame(sep_chn)) 

### Beta diversity

# The model implementations are NOT optimised for speed and can take a long time to run. 

D = pairwise(Jaccard(),Y,Y,dims = 1)
BC = pairwise(BrayCurtis(),props,props,dims = 1)

# conduct multi dimensional scaling and create output csv with mds coordinates
jcoords = nmds(D,2,80,0.0002,100000)[1]
jcoords ./= maximum(abs.(jcoords))
Gdata.jac1 =jcoords[:,1] 
Gdata.jac2 =jcoords[:,2]

bccoords = nmds(BC,2,80,0.0002,100000)[1]
bccoords ./= maximum(abs.(bccoords))
Gdata.bc1 =bccoords[:,1]
Gdata.bc2 =bccoords[:,2] 

CSV.write("denoised_mds.csv",Gdata)

# helper functions to sort data for the (bin)onmial model: creating matrices of (uni)ons and (int)ersections of the various communities and a
# "trt" matrix (indicating whether the respecitive communities are from the same nutrient treatment. 
bin_int(x,y) = sum(.*(x,y))
bin_uni(x,y) = sum(max.(x,y))

function bin_int(X)
    n= size(X,1)
    Y = Array{Int64}(undef,n,n)
    for i in 1:n
        for j in 1:n
            Y[i,j] = Y[j,i] = bin_int(X[i,:],X[j,:])
        end
    end
    Y
end
function bin_uni(X)
    n= size(X,1)
    Y = Array{Int64}(undef,n,n)
    for i in 1:n
        for j in 1:n
            Y[i,j] = Y[j,i] = bin_uni(X[i,:],X[j,:])
        end
    end
    Y
end

function bin_trt(X)
    n= size(X,1)
    Y = Array{Int64}(undef,n,n)
    for i in 1:n
        for j in 1:n
            Y[i,j] = Y[j,i] = (X[i] == X[j])
        end
    end
    Y
end

function get_binomial_data(Y,transfer,carbon)
    n_mask = transfer .== "none"
    l_mask = transfer .== "liquid"
    b_mask = transfer .== "bead"

    ((bin_uni(Y[n_mask,:]),bin_int(Y[n_mask,:]),bin_trt(carbon[n_mask])),   
    (bin_uni(Y[l_mask,:]),bin_int(Y[l_mask,:]),bin_trt(carbon[l_mask])),
    (bin_uni(Y[b_mask,:]),bin_int(Y[b_mask,:]),bin_trt(carbon[b_mask])))
end

binomial_data = get_binomial_data(Y,Gdata.transfer, Gdata.carbon)

# binomial model

@model function binmod(X,Y,trt,n = size(Y,1))
    a ~ Normal(0,4)
    b ~ Normal(0,4)
    λ ~ filldist(Normal(0,4),n)

    for j in 2:n
        for i in 1:j-1
            p = a*(1-trt[i,j]) + b * trt[i,j] + λ[i] + λ[j]
            Y[i,j] ~ BinomialLogit(X[i,j],p)
        end
    end
end

# wrapper model to perform inference on the 3 spatial structure treatments simultaneously.
@model function binmods(xn,xl,xb)
    @submodel prefix="sub1" N = binmod(xn[1],xn[2],xn[3])   
    @submodel prefix="sub2" l = binmod(xl[1],xl[2],xl[3])
    @submodel prefix="sub3" b = binmod(xb[1],xb[2],xb[3])
end


@model function normmod(Y,trt, trt2 = 1 .- trt)
    n = size(Y,1)
    a ~ Normal(0.5,1)
    b ~ Normal(0.5,1)
    σ ~ Exponential()
    λ ~ filldist(Normal(0.5,1),n)

    for j in 2:n
        for i in 1:j-1
            μ = a*trt2[i,j] + b * trt[i,j] + λ[i] + λ[j]
            Y[i,j] ~ Normal(μ,σ)
        end
    end
end
function get_normal_data(Y,transfer,carbon)
    n_mask = transfer .== "none"
    l_mask = transfer .== "liquid"
    b_mask = transfer .== "bead"
    ((Y[n_mask,n_mask],bin_trt(carbon)[n_mask,n_mask]),   
    (Y[l_mask,l_mask],bin_trt(carbon)[l_mask,l_mask]),
    (Y[b_mask,b_mask],bin_trt(carbon)[b_mask,b_mask]))
end

@model function normmods(xn,xl,xb)
    @submodel prefix="sub1" N = normmod(xn[1],xn[2])   
    @submodel prefix="sub2" l = normmod(xl[1],xl[2])
    @submodel prefix="sub3" b = normmod(xb[1],xb[2])
end

normal_data = get_normal_data(BC,Gdata.transfer, Gdata.carbon)

chn  = sample(binmods(binomial_data[1],binomial_data[2],binomial_data[3]),NUTS(0.65),MCMCThreads(),n_iter,n_chn,drop_warmup = true)[n_iter-n_keep+1:end]
chn_df = DataFrame(chn)

n_λ = vec(2mean(hcat([vec(chn["sub1.λ[$j]"]) for j in 1:36]...),dims=2))
l_λ = vec(2mean(hcat([vec(chn["sub2.λ[$j]"]) for j in 1:18]...),dims=2))
b_λ = vec(2mean(hcat([vec(chn["sub3.λ[$j]"]) for j in 1:17]...),dims=2))

chn_df.n_within = 1 .- Turing.logistic.(n_λ .+ vec(chn["sub1.b"]))
chn_df.n_without =1 .- Turing.logistic.(n_λ .+ vec(chn["sub1.a"]))
chn_df.l_within =1 .- Turing.logistic.(l_λ .+ vec(chn["sub2.b"]))
chn_df.l_without =1 .- Turing.logistic.(l_λ .+ vec(chn["sub2.a"]))
chn_df.b_within =1 .-Turing.logistic.(b_λ .+ vec(chn["sub3.b"]))
chn_df.b_without =1 .- Turing.logistic.(b_λ .+ vec(chn["sub3.a"]))
CSV.write("jaccard_chn.csv",chn_df) 


chn = sample(normmods(normal_data[1],normal_data[2],normal_data[3]),NUTS(0.65),MCMCThreads(),n_iter,n_chn,drop_warmup = true)[n_iter-n_keep+1:end]
chn_df = DataFrame(chn)

n_λ = vec(2mean(hcat([vec(chn["sub1.λ[$j]"]) for j in 1:36]...),dims=2))
l_λ = vec(2mean(hcat([vec(chn["sub2.λ[$j]"]) for j in 1:18]...),dims=2))
b_λ = vec(2mean(hcat([vec(chn["sub3.λ[$j]"]) for j in 1:17]...),dims=2))

chn_df.n_within =n_λ .+ vec(chn["sub1.b"])
chn_df.n_without =n_λ .+ vec(chn["sub1.a"])
chn_df.l_within =l_λ .+ vec(chn["sub2.b"])
chn_df.l_without =l_λ .+ vec(chn["sub2.a"])
chn_df.b_within =b_λ .+ vec(chn["sub3.b"])
chn_df.b_without =b_λ .+ vec(chn["sub3.a"])

CSV.write("bray_curtis_chn.csv",chn_df) 

beta_df = DataFrame([βpart(Y[Gdata.treatment .== trt,:]) for trt in unique(Gdata.treatment)])
bcdisp = dispersion(BC,Gdata.treatment)
beta_df[!,"Bray Curtis"] = round.(bcdisp.means, digits = 2)
beta_df = beta_df[[3,6,9,2,5,8,1,4,7],:]

trt_df = DataFrame()
trt_df[!,"Beads"] = vcat(fill("None",3),fill("New",3),fill("Same",3))
trt_df[!,"Nutrient"] = repeat(["Arabinose","Glucose","Leucine"],3)
beta_df = hcat(trt_df,beta_df)
CSV.write("beta_df.csv",beta_df)

