#load packages
using TikzPictures, FileIO
using CSV, DataFrames, Turing,DataFrameMacros,
 Colors,Random,CairoMakie
 CairoMakie.activate!()
 using Rsvg
 using Cairo
 Random.seed!(123)
 n_gen = 17 # total number of genera
 include("DiversityFunctions.jl")
 #Use the colour palette from the manuscript.
 nobeadclr = RGB(([0,47,74] ./255)...)
 liquiddclr = RGB(([214,42,41] ./255)...)
 beadclr = RGB(([253,192,76] ./255)...)
 palette = [beadclr,liquiddclr,nobeadclr]
Gdata =CSV.read("denoised_mds.csv", DataFrame)

X = Array(Gdata[:,2:n_gen]) # Array of genus counts
Gdata.pielou = pielou.(eachrow(X))
Y = X .> 0 # Presence/absence matrix
props = X ./ Gdata.depth # proportions
n_row,n_genera = size(X)
n_treat = length(unique(Gdata.treatment))
source_inds = [findfirst(x -> x == t,["ara","glu","leu"]) for t in Gdata.carbon];
transfer_inds = [findfirst(x -> x == t,["none","liquid","bead"]) for t in Gdata.transfer];
### alpha diversity

n_iter = 2000 # number of MCMC samples to take per chain
n_keep = 1000 # number of MCMC samples to keep per chain
n_chn = 4 # MCMC chains per model
page_width = 6.5 # inches

size_in_inches = (page_width, 4)
dpi = 2700
size_in_pixels = size_in_inches .* dpi

Random.seed!(123)

function plot_lab(fig)
    axs = fig.content
    filter!(x -> typeof(x) == Axis,axs)
    n = length(axs)
    letters = collect('a':'z')[1:n]
  for (label, ax) in zip(letters, axs)
      ax.subtitle = "("*label*")"
  end
end
word_font = "Helvetica"
naked_theme = Theme(
    Axis = (
        titlealign = :left,
        yticklabelsvisible = false,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = true,
        topspinevisible = false,
        xgridvisible   = false,
        ygridvisible   = false,
        yticksvisible   = false,
        xlabelsize = 20,
        titlefont = word_font,
        xticklabelfont = word_font,
        yticklabelfont = word_font,
        ylabelfont = word_font,
        xlabelfont = word_font,
        ylabelsize = 20,
    ),
    Legend = (
      titlefont = word_font,
      labelfont = word_font
    )
)

basic_theme = Theme(
    Axis = (
        titlealign = :left,
        xgridvisible   = false,
        ygridvisible   = false,
        xlabelsize = 20,
        titlefont = word_font,
        xticklabelfont = word_font,
        yticklabelfont = word_font,
        ylabelfont = word_font,
        xlabelfont = word_font,
        ylabelsize = 20,
    ),
    Legend = (
      titlefont = word_font,
      labelfont = word_font
    )
)

set_theme!(naked_theme);
#### Figure 1: model

tp = include("model_scheme.jl")

TikzPictures.save(SVG("model_scheme"), tp)

filename_in = "model_scheme.svg"
filename_out = "model_scheme.png"

Rsvg.set_default_dpi(300.0)
r = Rsvg.handle_new_from_file(filename_in);
d = Rsvg.handle_get_dimensions(r);
cs = Cairo.CairoImageSurface(d.width,d.height,Cairo.FORMAT_ARGB32);
c = Cairo.CairoContext(cs);
Rsvg.handle_render_cairo(c,r);
Cairo.write_to_png(cs,filename_out);

img = load("model_scheme.png")
fig = Figure()
ax1 = Axis(fig[1:2,1:2])
CairoMakie.image!(ax1,rotr90(img))
ax1.aspect = DataAspect()
env = CSV.read("model_res_2.csv",DataFrame)

df =  @subset(env, :richness > 0)
static_df = @subset(df, :dispersal ==0, :c ==0.9)
shaken_df = @subset(df, :dispersal >0, :c ==0.9)

ax2 = Axis(fig[1,3])
bead = hist!(ax2,static_df.richness,bins = 1:100, strokewidth = 1,strokecolor = (:black, 0.5), color = (beadclr,1.0), normalization = :pdf)
liq = hist!(ax2,shaken_df.richness,bins = 1:100,strokewidth = 1,strokecolor = (:black, 0.5), color = (nobeadclr,0.7),normalization = :pdf)

# for beta diversity metrics use communities with atleast 2 species

df =  @subset(df, :richness > 1)
static_df = @subset(df, :dispersal ==0, :c ==0.9)
shaken_df = @subset(df, :dispersal >0, :c ==0.9)

ax3 = Axis(fig[2,3])
bead = hist!(ax3,static_df.pie,bins = 0.0:0.05:1, strokewidth = 1,strokecolor = (:black, 0.5), color = (beadclr,1.0), normalization = :pdf)
liq = hist!(ax3,shaken_df.pie,bins = 0.0:0.05:1,strokewidth = 1,strokecolor = (:black, 0.5), color = (nobeadclr,0.7),normalization = :pdf)

ax4 = Axis(fig[3,1])
bead = hist!(ax4,static_df.jac,bins = 0.0:0.05:1, strokewidth = 1,strokecolor = (:black, 0.5), color = (beadclr,1.0), normalization = :pdf)
liq = hist!(ax4,shaken_df.jac,bins = 0.0:0.05:1,strokewidth = 1,strokecolor = (:black, 0.5), color = (nobeadclr,0.7),normalization = :pdf)

ax5 = Axis(fig[3,2])
bead = hist!(ax5,static_df.envdet,bins = 0.0:0.05:1, strokewidth = 1,strokecolor = (:black, 0.5), color = Makie.to_color((beadclr,1.0)), normalization = :pdf)
liq = hist!(ax5,shaken_df.envdet,bins = 0.0:0.05:1,strokewidth = 1,strokecolor = (:black, 0.5), color = Makie.to_color((nobeadclr,0.7)),normalization = :pdf)
Legend(fig[3, 3], [bead, liq], ["Spatial", "Well mixed"], linecolor = :red,tellheight = false, tellwidth = false)

xlims!(ax2,(0,100))
xlims!(ax3,(0,1))
xlims!(ax4,(0,1))
xlims!(ax5,(0,1))

ax2.xlabel = "Species richness"
ax3.xlabel = "Pielou's Index"
ax4.xlabel = "Jaccard's Index"
ax5.xlabel = "Jaccard's Index"
ax2.aspect = 2
ax3.aspect = 2
ax4.aspect = 2
ax5.aspect = 2
ax1.bottomspinevisible = false
hidedecorations!(ax1)
plot_lab(fig)
colgap!(fig.layout, 28)
rowgap!(fig.layout, 2)

CairoMakie.save("model_res.png",fig, px_per_unit = 3);


set_theme!(basic_theme);

#### Figure 2: proportions

## helper funcs
function poly3(t, p0, p1, p2, p3)
    Point2f((1-t)^3 .* p0 .+ t*p1*(3*(1-t)^2) + p2*(3*(1-t)*t^2) .+ p3*t^3)
end
function BezierPath(o, f, co, cf; t = range(0,1, length=30))
    return [poly3(t, o, co, cf, f) for t in t]
end

function supLine(p1, p2; x=0,y=8)
    [p1 .+ Point2f(x,y), p1, p2, p2 .+ Point2f(x,y)]
end
function supLine2(p1, p2; x=0,y=8)
    [p1 .- Point2f(x,y), p1, p2, p2 .- Point2f(x,y)]
end

## plot
cmap = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", 
"#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", 
"#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000"]

f = Figure(; resolution = (1200,800))
ax = Axis(f[1,1:9], ylabel = "Proportion of reads", xlabel = "")

hidedecorations!(ax)
ylims!(ax,-0.3,1.1)
xlims!(ax,0,72)
hidespines!.(ax)
    
    p = vcat([props[Gdata.treatment .== trt,:] for trt in ["ara_none","glu_none", "leu_none",
                                                        "ara_liquid","glu_liquid", "leu_liquid",
                                                        "ara_bead","glu_bead", "leu_bead"]]...)
    tbl = (x = vcat([fill(i,size(p,2)) for i in 1:size(p,1)]...),
    height = vcat(eachrow(p)...),
    grp = vcat([collect(1:size(p,2)) for i in 1:size(p,1)]...),
    )
    barplot!(ax,tbl.x, tbl.height,
            stack = tbl.grp,
            color = tbl.grp,
            offset = 15,
            colormap = cmap[1:size(p,2)]
            )

ops = [[[18,45,63.5][i],-0.05] for i in 1:3]
fps = [[[30,35.5,41][i],-0.2] for i in 1:3]

    supls = [supLine(Point2f([[1,37,55][i],-0.05]),
    Point2f([[36,54,71][i],-0.05]),y = 0.03) for i in 1:3]
[lines!(ax, supls[k], color = reverse(palette)[k]) for k in 1:3]
[lines!(ax, BezierPath(ops[k], fps[k], [ops[k][1],ops[k][2]-0.1],
    [fps[k][1],fps[k][2]+0.1]); color = reverse(palette)[k]) for k in 1:3]
text!(ax,[30,35.5,41], [-0.24,-0.24,-0.24], text = ["No\nbeads", "New\nbeads", "Same\nbeads"],
    font = word_font, align = (:center, :center))

text!([6.5,18.5,30.5,39.5,45.5,51.5,57,62.5,68.5],
        fill(1.08,9), text = repeat(["ara","glu","leu"],3), 
        font = word_font,align = (:center, :center))

        supls = [supLine2(Point2f([[1.0,13,25,37,43,49,55,60,66][i],1.05]),
                Point2f([[12.0,24,36,42,48,54,59,65,71][i],1.05]),y = 0.03) for i in 1:9]
                
[lines!(ax, supls[k], color = :black) for k in 1:9]

Legend(f[1, 10],
[MarkerElement(color = cmap[c], marker = :rect, markersize = 15) for c in 1:size(p,2)],
names(Gdata)[2:n_gen])

f.attributes.resolution = size_in_pixels
CairoMakie.save("props.png",f, px_per_unit = 3)

#### Figure 3: alpha

fig = Figure()
box = Axis(fig[1,1:2])
box.ylabel = "Genus richness"
cx = [findfirst(x -> x == c, unique(Gdata.carbon)) for c in Gdata.carbon]
tx = [findfirst(x -> x == t, reverse(unique(Gdata.transfer))) for t in Gdata.transfer]
x = tx .+ 5cx 
jitter = rand(truncated(Normal(0,0.2),-0.5,0.5),length(x))
yjitter = rand(truncated(Normal(0,0.1),-0.1,0.1),length(x))

boxplot!(box,x, Gdata.richness, color = reverse(palette)[tx],strokewidth = 2, show_outliers = false)
scatter!(box,x .+jitter, Gdata.richness .+ yjitter, color = (:white,0.8),strokewidth = 2,strokecolor = reverse(palette)[tx])

chn_df =CSV.read("richness_chn.csv",DataFrame) 
chn2 = Chains(Array(chn_df), names(chn_df))

# plot chains
chn_datar = Array(DataFrame(group(chn2,:muy)))[:,3:end]
alphadiv = Axis(fig[1,3])
alphadiv.ylabel = "Genus richness"

n =violin!(alphadiv,fill(1,n_keep * n_chn),chn_datar[:,1],color =(nobeadclr,0.7),strokewidth = 2)
l = violin!(alphadiv,fill(2,n_keep * n_chn),chn_datar[:,2],color =(liquiddclr,0.7),strokewidth = 2)
b = violin!(alphadiv,fill(3,n_keep * n_chn),chn_datar[:,3],color =(beadclr,0.7),strokewidth = 2)
alphadiv.yaxisposition = :right
linkyaxes!(box,alphadiv)
hidexdecorations!.(fig.content)

alphadiv.yticks = 8:17
box.yticks = 8:17

chn_df =CSV.read("pielou_chn.csv",DataFrame) 
sep_chn = Chains(Array(chn_df), names(chn_df))


box = Axis(fig[2,1:2])
box.ylabel = "Pielou's evenness"

cx = [findfirst(x -> x == c, unique(Gdata.carbon)) for c in Gdata.carbon]
tx = [findfirst(x -> x == t, reverse(unique(Gdata.transfer))) for t in Gdata.transfer]
x = tx .+ 5cx 
jitter = rand(truncated(Normal(0,0.02),-0.05,0.05),length(x))
yjitter = rand(truncated(Normal(0,0.01),-0.01,0.01),length(x))

boxplot!(box,x, Gdata.pielou, color = reverse(palette)[tx],strokewidth = 2, show_outliers = false)
scatter!(box,x .+jitter, Gdata.pielou .+ yjitter, color = (:white,0.8),strokewidth = 2,strokecolor = reverse(palette)[tx])


chn_datar = hcat(vcat(sep_chn["p[1]"]...),vcat(sep_chn["p[2]"]...),vcat(sep_chn["p[3]"]...))

# plot chains
alphadiv = Axis(fig[2,3], xlabel = "μ")
n =violin!(alphadiv,fill(1,n_keep * n_chn),chn_datar[:,1],color =Makie.to_color((nobeadclr,0.7)),strokewidth = 2)
l = violin!(alphadiv,fill(2,n_keep * n_chn),chn_datar[:,2],color =Makie.to_color((liquiddclr,0.7)),strokewidth = 2)
b = violin!(alphadiv,fill(3,n_keep * n_chn),chn_datar[:,3],color =Makie.to_color((beadclr,0.7)),strokewidth = 2)
alphadiv.yaxisposition = :right
alphadiv.ylabel = "Pielou's evenness"
linkyaxes!(box,alphadiv)
hideydecorations!.(fig.content, label = false,ticks = false,ticklabels = false)
hidexdecorations!.(fig.content, ticks = false,ticklabels = false)

box.xticks =([7,12,17],["Arabinose", "Glucose", "Leucine"])
alphadiv.xticks =([1,2,3],["No \n beads", "New \n beads", "Same \n beads"])

plot_lab(fig)
alphadiv.xticklabelsize =20
box.xticklabelsize =20

axislegend(box, [n,l,b], ["No beads", "New beads", "Same beads"], position = :rt,framevisible = false)
CairoMakie.save("alpha.png",fig, px_per_unit = 3)


#### Figure 4: beta
function get_centroids(coords,df, transfer)
    sub_coords = coords[df.transfer .== transfer, :]
    sub_carbon = df.carbon[df.transfer .== transfer]
    centre = Point(Tuple(mean(sub_coords, dims = 1)))
    ara = Point(Tuple(mean(sub_coords[sub_carbon .== "ara",:], dims = 1)))
    glu = Point(Tuple(mean(sub_coords[sub_carbon .== "glu",:], dims = 1)))
    leu = Point(Tuple(mean(sub_coords[sub_carbon .== "leu",:], dims = 1)))
    return centre, ara, glu, leu
end
jcoords= hcat(Gdata.jac1,Gdata.jac2)
bccoords= hcat(Gdata.bc1,Gdata.bc2)

figy = Figure()
cmap = [(p,0.7) for p in palette]
markers =[[:utriangle,:circle,:diamond][i] for i in source_inds]
jac = Axis(figy[1,2])
scatter!(jac,Gdata.jac1, Gdata.jac2; 
    color =transfer_inds , colormap = reverse(cmap), marker = markers)
    bead_centroids = get_centroids(jcoords, Gdata, "bead")
    liquid_centroids = get_centroids(jcoords, Gdata, "liquid")
    none_centroids = get_centroids(jcoords, Gdata, "none")
    [lines!(jac,[bead_centroids[[1,i]]...],linestyle = :dash, color = beadclr) for i in 2:4]
    [lines!(jac,[liquid_centroids[[1,i]]...],linestyle = :dash, color = liquiddclr) for i in 2:4]
    [lines!(jac,[none_centroids[[1,i]]...],linestyle = :dash, color = nobeadclr) for i in 2:4]

    [scatter!(jac,bead_centroids[i], color = beadclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
    [scatter!(jac,liquid_centroids[i], color = liquiddclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
    [scatter!(jac,none_centroids[i], color = nobeadclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
    

jac.aspect=DataAspect()
mn,mx = minimum(jcoords)-0.2, maximum(jcoords)+0.2
jac.limits= (mn,mx,mn,mx)
jac.ylabel = "Dim 2"
hidexdecorations!(jac)

bc = Axis(figy[2,2])
scatter!(bc,Gdata.bc1, Gdata.bc2;
    color =transfer_inds , colormap = reverse(cmap), marker = markers)
    bead_centroids = get_centroids(bccoords, Gdata, "bead")
    liquid_centroids = get_centroids(bccoords, Gdata, "liquid")
    none_centroids = get_centroids(bccoords, Gdata, "none")
    [lines!(bc,[bead_centroids[[1,i]]...],linestyle = :dash, color = beadclr) for i in 2:4]
    [lines!(bc,[liquid_centroids[[1,i]]...],linestyle = :dash, color = liquiddclr) for i in 2:4]
    [lines!(bc,[none_centroids[[1,i]]...],linestyle = :dash, color = nobeadclr) for i in 2:4]
    [scatter!(bc,bead_centroids[i], color = beadclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
    [scatter!(bc,liquid_centroids[i], color = liquiddclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
    [scatter!(bc,none_centroids[i], color = nobeadclr,strokewidth = 1, marker = [:rect,:utriangle,:circle,:diamond][i]) for i in 1:4] 
bc.xlabel = "Dim 1"
bc.ylabel = "Dim 2"

bc.aspect=DataAspect()
mn,mx = minimum(bccoords)-0.2, maximum(bccoords)+0.2
bc.limits= (mn,mx,mn,mx)

markers = [:rect,:utriangle,:circle,:diamond]

mk_els =[MarkerElement(color = :black, marker = mk, markersize = 15) for mk in markers]
clr_els = [MarkerElement(color = clr, marker = :rect, markersize = 15) for clr in reverse(palette)]
lin_el = LineElement(color = :black, linestyle = :dash)

    figy[1:2,1] =Legend(figy,
    [vcat(mk_els...,lin_el), clr_els],
    [["Culture method centroid","Arabinose","Glucose","Leucine","Distance to centroid"], ["No beads", "New beads", "Same beads"]],
    ["Marker", "Colour"])


chn_df =CSV.read("jaccard_chn.csv",DataFrame) 
chn = Chains(Array(chn_df), names(chn_df))


n_λ = vec(2mean(hcat([vec(chn["sub1.λ[$j]"]) for j in 1:36]...),dims=2))
l_λ = vec(2mean(hcat([vec(chn["sub2.λ[$j]"]) for j in 1:18]...),dims=2))
b_λ = vec(2mean(hcat([vec(chn["sub3.λ[$j]"]) for j in 1:17]...),dims=2))

ax3 = Axis(figy[1,3])
violin!(ax3,fill(1,n_keep * n_chn),1 .- Turing.logistic.(n_λ .+ vec(chn["sub1.b"])), color = (nobeadclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(1,n_keep * n_chn),1 .- Turing.logistic.(n_λ .+ vec(chn["sub1.a"])), color = (nobeadclr,0.9),side = :right, strokewidth = 2)

violin!(ax3,fill(1.5,n_keep * n_chn),1 .- Turing.logistic.(l_λ .+ vec(chn["sub2.b"])), color = (liquiddclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(1.5,n_keep * n_chn),1 .- Turing.logistic.(l_λ .+ vec(chn["sub2.a"])), color = (liquiddclr,0.9),side = :right, strokewidth = 2)

violin!(ax3,fill(2,n_keep * n_chn),1 .-Turing.logistic.(b_λ .+ vec(chn["sub3.b"])), color = (beadclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(2,n_keep * n_chn),1 .- Turing.logistic.(b_λ .+ vec(chn["sub3.a"])), color = (beadclr,0.9),side = :right, strokewidth = 2)
ax3.aspect = 1
hidexdecorations!(ax3)
ax3.ylabel = "Jaccard's dissimilarity"


chn_df =CSV.read("bray_curtis_chn.csv",DataFrame) 
chn = Chains(Array(chn_df), names(chn_df))

n_λ = vec(2mean(hcat([vec(chn["sub1.λ[$j]"]) for j in 1:36]...),dims=2))
l_λ = vec(2mean(hcat([vec(chn["sub2.λ[$j]"]) for j in 1:18]...),dims=2))
b_λ = vec(2mean(hcat([vec(chn["sub3.λ[$j]"]) for j in 1:17]...),dims=2))

ax3 = Axis(figy[2,3])
violin!(ax3,fill(1,n_keep * n_chn),n_λ .+ vec(chn["sub1.b"]), color = (nobeadclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(1,n_keep * n_chn),n_λ .+ vec(chn["sub1.a"]), color = (nobeadclr,0.9),side = :right, strokewidth = 2)

violin!(ax3,fill(1.5,n_keep * n_chn),l_λ .+ vec(chn["sub2.b"]), color = (liquiddclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(1.5,n_keep * n_chn),l_λ .+ vec(chn["sub2.a"]), color = (liquiddclr,0.9),side = :right, strokewidth = 2)

violin!(ax3,fill(2,n_keep * n_chn),b_λ .+ vec(chn["sub3.b"]), color = (beadclr,0.9),side = :left, strokewidth = 2)
violin!(ax3,fill(2,n_keep * n_chn),b_λ .+ vec(chn["sub3.a"]), color = (beadclr,0.9),side = :right, strokewidth = 2)

plot_lab(figy)
ax3.aspect = 1
ax3.xticks = ([1,1.5,2], ["No beads", "New beads", "Same beads"])
ax3.xticklabelrotation = pi/3
ax3.xticklabelsize = 20
ax3.ylabel = "Bray-Curtis dissimilarity"
figy.attributes.resolution = size_in_pixels
CairoMakie.save("beta.png",figy,px_per_unit = 3)

#### suppl chem

tp =TikzPicture(L"""
\newsavebox\arabinose
\sbox\arabinose{
\chemfig[bond style={line width=2pt}, double bond sep = 3pt]{ 
           HO% 7
     -[,,2]% 5
    -[:300]% 4
              (
    -[:240,,,2]HO% 8
              )
          -% 3
              (
    -[:300,,,1]OH% 9
              )
     -[:60]% 2
              (
        -[,,,1]OH% 10
              )
    -[:120]% 1
    -[:180]O% 6
              (
        -[:240]% -> 5
              )
                  }
}
\newsavebox\glucose
\sbox\glucose{
\chemfig[bond style={line width=2pt}, double bond sep = 3pt]{ 
               HO% 8
      -[:60,,2]% 6
              -% 5
                  (
        -[:300,,,1]OH% 9
                  )
         -[:60]% 4
                  (
            -[,,,1]OH% 10
                  )
        -[:120]% 3
                  (
         -[:60,,,1]OH% 11
                  )
        -[:180]% 2
                  (
            -[:240]O% 7
            -[:300]% -> 6
                  )
        -[:120]% 1
    -[:180,,,2]HO% 12
    }
}
\newsavebox\leucine
\sbox\leucine{
\chemfig[bond style={line width=2pt}, double bond sep = 3pt]{               % 1
         -[:90]% 2
                  (
            -[:150]% 3
                  )
         -[:30]% 4
        -[:330]% 5
                  (
        -[:270,,,1]NH_2% 9
                  )
         -[:30]% 6
                  (
             =[:90]O% 7
                  )
    -[:330,,,1]OH% 8
    }
}
\node[text=black,inner sep=0pt, outer sep=0pt, scale=  1.12,label =below:{Arabinose}] at (0.8,0.1) (A) {};
\node[text=black,inner sep=0pt, outer sep=0pt, scale=  1.12,label =above:{Glucose}] at (2.3,3.1) (G) {};
\node[text=black,inner sep=0pt, outer sep=0pt, scale=  1.12,label = below:{Leucine}] at (3.8,0.1) (L) {};
\node [below =0.45cm of L, scale = 0.2] {};
\node [above=0.4cm of G, scale = 0.2] {};
\node [below =0.45cm of A, scale = 0.2] {};
\draw [thick]  (A) -- node [above, sloped] {$0.68$} (G);
\draw [thick]  (A) -- node [above, sloped] {$0.1$} (L);
\draw [thick]  (G) -- node [above, sloped] {$0.08$} (L);
""", preamble = "\\usepackage{amsmath}
\\definecolor{nobeadclr}{RGB}{0,47,74}
\\definecolor{liquidclr}{RGB}{214,42,41}
\\definecolor{beadclr}{RGB}{253,192,76}
\\definecolor{Ahair}{RGB}{244,237,70}
\\definecolor{Adress1}{RGB}{180,242,255}
\\definecolor{Adress2}{RGB}{92,225,244}
\\definecolor{Askin}{RGB}{255,234,202}
\\definecolor{Rcollar}{RGB}{255,65,65}
\\usepackage{chemfig}
\\usepackage[scaled]{helvet}
\\usepackage[T1]{fontenc}
\\renewcommand\\familydefault{\\sfdefault}
\\usepackage{pgfplots}
\\usepackage{tikz-network}
"
)


TikzPictures.save(SVG("chem"), tp)

filename_in = "chem.svg"
filename_out = "chem.png"

Rsvg.set_default_dpi(300.0)
r = Rsvg.handle_new_from_file(filename_in);
d = Rsvg.handle_get_dimensions(r);
cs = Cairo.CairoImageSurface(d.width,d.height,Cairo.FORMAT_ARGB32);
c = Cairo.CairoContext(cs);
Rsvg.handle_render_cairo(c,r);
Cairo.write_to_png(cs,filename_out);

