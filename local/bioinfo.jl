# load packages
using DataFrames, CSV


reorder =  sortperm(parse.(Int,sort(string.(1:132))))

liq_inds = 85:102
bead_inds = 103:120
tsb_inds = 121:132
none_inds = 1:84
sing_inds = 1:36
multi_inds = 37:84
sample_names = collect(1:132)
seqtab = CSV.read("seq_all_q.csv",DataFrame)[reorder,2:end]# 1st column just contains row numbers
taxa = CSV.read("taxa_all_q.csv",DataFrame)
df = DataFrame([sample_names,seqtab[:,1]], ["sample",taxa[1,6]] )
for (i,genus) in enumerate(taxa[2:end,7])
    if (genus != "NA")
        
        if in(genus,names(df))
            df[!,genus] .+= seqtab[:,i]
        else
            df[!,genus] =  seqtab[:,i]
            println(genus)
        end
    else
    end
end

n_gen = size(df)[2]

df[!,"total"] =vec(sum(Array(df[:,2:end]),dims = 2))
df[!,"transfer"] .= "none"
df.transfer[bead_inds] .= "bead"
df.transfer[liq_inds] .= "liquid"
df.transfer[liq_inds] .= "liquid"
df = df[vcat(sing_inds,liq_inds,bead_inds),:]
df.carbon = vcat([fill("ara",12), fill("glu",12), fill("leu",12),
            fill("ara",6), fill("glu",6), fill("leu",6),
            fill("ara",6), fill("glu",6), fill("leu",6)]...)

df[!,"treatment"] = df.carbon .* "_" .* df.transfer

Gdata = sort(df,:treatment)

for i in 2:n_gen
    if sum(Gdata[!,i]) < 10
     Gdata[!,i] .= 0
     n_gen -=1
    end
 end

 Gdata = Gdata[!, .!all.(==(0), eachcol(Gdata))]
    Gdata.richness = sum(eachcol(Gdata[:,2:n_gen] .>0))
    Gdata.depth = sum(eachcol(Gdata[:,2:n_gen]))
    # remove single outlier with very low read depth
    Gdata = Gdata[Gdata.depth .> minimum(Gdata.depth),:]
    max_depth = maximum(Gdata.depth)

CSV.write("denoised_q.csv",Gdata)



