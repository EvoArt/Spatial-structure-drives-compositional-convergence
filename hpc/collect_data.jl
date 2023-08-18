using Pkg
Pkg.activate(".")
using CSV, DataFrames

data = filter(x -> contains(x,res"[0-9]+"),readdir())
data = dropmissing(CSV.read(data, DataFrame))
CSV.write("model_res.csv",data)

data = filter(x -> contains(x,r"specmat[0-9]+"),readdir())
data = dropmissing(CSV.read(data, DataFrame))
CSV.write("model_taxa.csv",data)
