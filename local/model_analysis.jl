using CSV,DataFrames,StatsBase,DataFrameMacros, 
Distances, BetaDispersion
include("DiversityFunctions.jl")

env = CSV.read("model_res.csv",DataFrame)
env.richness .-= 1
specenv = Array{Int}(CSV.read("model_taxa.csv",DataFrame))
specenv2 = specenv .>0

jac = vec(vcat([dispersion(specenv2[1+(i-1)*10:i*10,:],ones(10), Jaccard).residuals[1] for i in 1:size(env)[1] ÷ 10]...))
envdet = vcat([repeat(vcat([evaluate(Jaccard(),specenv2[j,:],specenv2[j+10,:]) for j in i:i+9]...),2) for i in 1:20:size(env)[1]]...)
env.jac = vec(jac)
env.envdet = vec(envdet)

env.pie = pielou.(eachrow(Array{Int}(specenv)))
env.dispersal = env.dispersal .>0

CSV.write("model_res_2.csv",env)

